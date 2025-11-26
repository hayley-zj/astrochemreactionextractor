import requests
import re
import json
import time
import configparser
from astrochemreactionextractor import logger
import pandas as pd

from astrochemreactionextractor.utils.utils_reaction import read_and_split_txt, add_rate_coefficient


config = configparser.ConfigParser()
# 读取配置文件
config.read('./astrochemreactionextractor/config.ini')
# 读取字符串值
llm_url = config.get('api', 'url')
llm_api_key = config.get('api', 'api_key')

def has_chinese(text):
    return bool(re.search('[\u4e00-\u9fff]', text))


def make_post_request(url, headers, payload, max_retries=3, retry_delay=1):
    retries = 0
    while retries < max_retries:
        try:
            response = requests.post(
                url, headers=headers, json=payload, timeout=600)
            return response
        except:
            print(f"{url} Retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)
        retries += 1
    logger.error(f"Interface request failed: {url}")
    return None


def deepseek_api_reactions_pdf(txt_contents):
    api_key = llm_api_key
    url = llm_url
    headers = {
        'Authorization': f'Bearer {api_key}',
        'Content-Type': 'application/json'
    }
    temperature = 0.01
    max_tokens = 12000

    txt_reactions = []
    prompt_tokens = 0
    completion_tokens = 0
    formula_dict = {}

    for index, row in txt_contents.iterrows():
        reactions_txt = row['reactions_txt']
        col2_coordinates = row['coordinates']
        text_chunks = read_and_split_txt(reactions_txt, max_tokens)
        start_text = ''
        chunk_temp = ''
        rpos_table = 0
        pos_table = 0
        for chunk in text_chunks:
            if chunk.count("→") > 10 and "Table" in chunk:
                if chunk.count("Table ") > 1:
                    rpos_table = chunk.rfind("Table ")
                    start_text = chunk[rpos_table:rpos_table + 110] + '\n\n'
                else:
                    pos_table = chunk.find("Table ")
                    if pos_table > 200:
                        chunk = start_text + chunk + '\n\n'

                    start_text = chunk[len(
                        start_text) + pos_table:len(start_text) + pos_table + 110] + '\n\n'
            if start_text != '' and "Table" not in chunk[0:200]:
                chunk = start_text + chunk + '\n\n'
            payload = {
                "model": 'deepseek-v3',
                "max_tokens": max_tokens,
                "temperature": temperature,
                "messages": [
                    {
                        "role": "user",
                        "content":
                                f"""
                                As a chemistry expert, please accurately extract the following information from the text:

                                1. **Chemical reaction formula** (correct garbled characters)
                                2. **Parameters related to chemical reaction rate coefficients, such as α, β, γ, EA (or Ea), rate/reaction coefficients, rate constant, temperature, pressure, etc. **
                                3. **Calculation formula for chemical reaction rate coefficient**
                                4. **Reaction type and reaction location** (only if explicitly mentioned in the original text)

                                ### Extraction rules

                                1. Chemical reaction extraction:

                                - Splicing dispersed molecular formulas

                                - `/` splits multiple reactions (e.g. `A/B + C → D/E` → `A + C -> D` and `B + C -> E`)

                                - The space before `->` is considered as no reactant (e.g. `\n->X` → `-> X`)

                                - Only extract the reaction formulas directly given in the original text

                                2. Extraction of information related to chemical reaction rate coefficient:
                                
                                - Variable actual values: Only variables with actual values ​​are extracted.

                                - Reaction rate coefficient calculation formulas (including variable definitions). Only reaction rate coefficient calculation formulas are extracted; other types of formulas are not extracted.

                                - The corresponding object of the reaction rate coefficient calculation formula. For example, if the rate coefficient calculation formula corresponds to the reaction type, the extraction result is '1:k = αζ'. For example, if the rate coefficient calculation formula in 'The rate coefficient of reaction C3H9+ + e− -> H + C3H8 is given by k = α(T/300)^β e^(-γ/T)' corresponds to the chemical reaction formula, the corresponding field is output directly after the reaction formula.

                                - **Formula formatting**:

                                - Exponentiation with `^` (`10 7` → `10^7`)

                                - Use `_` for subscripts (`A V` → `A_v`)

                                - **Not Extracted**:

                                - Isolated values ​​with undefined formulas

                                - Ambiguous references ("see Tab", etc.)

                                3. Reaction type and reaction location:

                                - **Reaction Type**: Extract chemical reaction types such as 'dissociative recombination reactions';

                                - **Reaction position**: Common reaction positions include 'ice', 'dust', 'grain', 'gas', 'Surface', etc.;

                                ### Output format

                                ```plaintext

                                # Chemical reaction formula lines (one line per formula, missing fields are not output)

                                [Corrected reaction formula]; reaction type: [Reaction type]; reaction position: [Reaction position]; variable values: [Variable name: actual value]; calculation formula: [Calculation formula for reaction rate coefficient], [Variable definition]; 

                                # Rate formula line (one line per formula, this field is not required)
                                [Reaction Type]: [Reaction Rate Coefficient Calculation Formula], [Variable Definition]

                                ```

                                ### Strict ban

                                1. No speculation about what is not mentioned

                                2. Adding Chinese/explanatory text is prohibited

                                3. Disable output of example content

                                ### Please start processing the following text: {chunk}
                                """
                    }
                ]
            }

            response = make_post_request(url, headers, payload)

            if response != None and response.status_code == 200:
                data = json.loads(response.text)
                if 'choices' not in data:
                    logger.warning(f"No result in data, chunk: {chunk}")
                    continue
                result_text = data['choices'][0]['message']['content'].replace(
                    '。', '')
                prompt_tokens += data['usage']['prompt_tokens']
                completion_tokens += data['usage']['completion_tokens']
                reaction_list = [reaction.strip() for reaction in result_text.split(
                    '\n') if reaction.strip()]
                for recation in reaction_list:
                    reaction = recation.replace('⇌', '->')
                    reaction = re.sub(r'^\d+(?:\,|\.)\s*', '', reaction)
                    if not has_chinese(reaction):
                        if '->' in reaction or '→' in reaction or '⇌' in reaction:
                            recation_split = reaction.split(';', 1)
                            if len(recation_split) > 1:
                                if 'EA:' in recation_split[1] and 'reaction position' not in recation_split[1]:
                                    txt_reactions.append(
                                        [recation_split[0], recation_split[1] + ', reaction position: surface', col2_coordinates])
                                else:
                                    txt_reactions.append(
                                        [recation_split[0], recation_split[1], col2_coordinates])
                            else:
                                txt_reactions.append(
                                    [recation_split[0], '', col2_coordinates])
                        if 'calculation formula is:' in reaction:
                            try:
                                remaining_text = reaction.split(
                                    ':', 1)[1].strip()
                                parts = remaining_text.split(
                                    ':', 1)
                                formula_dict.update(
                                    {parts[0]: parts[1].strip()})
                            except:
                                continue
            else:
                logger.error(f"Request failed!!!")

    df_retrieval_txt_reactions = pd.DataFrame(
        txt_reactions, columns=['reactions_wx', 'params', 'coordinates'])

    df_results = add_rate_coefficient(df_retrieval_txt_reactions, formula_dict)
    return df_results, prompt_tokens, completion_tokens


def deepseek_api_formula(chemical_reaction):
    api_key = llm_api_key
    url = llm_url

    headers = {
        'Authorization': f'Bearer {api_key}',
        'Content-Type': 'application/json'
    }
    temperature = 0.1
    max_tokens = 12000

    results = {}
    payload = {
        "model": 'deepseek-v3',
        "max_tokens": max_tokens,
        "temperature": temperature,
        "messages": [
            {
                "role": "user",
                "content":
                        f"""
                        Task Description:
                        Your task is to extract all chemical molecular formulas from the chemical reaction formula. For example, the molecular formulas in 'CH + CH2CHCHCH2 -> c-C5H6 + H' are 'CH', 'CH2CHCHCH2', 'c-C5H6', and 'H'. Please process the chemical reaction formula and extract all the chemical molecular formulas from it:
                        {chemical_reaction}

                        Output requirements:
                        Output the result directly without any extra text, One line is a chemical formula.
                        """
            }
        ]
    }

    response = make_post_request(url, headers, payload)
    result_text = ''

    prompt_tokens = 0
    completion_tokens = 0
    if response != None and response.status_code == 200:
        data = json.loads(response.text)
        if 'choices' in data:
            result_text = data['choices'][0]['message']['content']
            prompt_tokens = data['usage']['prompt_tokens']
            completion_tokens = data['usage']['completion_tokens']
    else:
        logger.warn(f"Request failed!")

    return result_text, prompt_tokens, completion_tokens


def deepseek_api_species(chemical_reaction, txt_contents):
    api_key = llm_api_key
    url = llm_url
    headers = {
        'Authorization': f'Bearer {api_key}',
        'Content-Type': 'application/json'
    }
    temperature = 0.1
    max_tokens = 12000

    if txt_contents == '':
        txt_contents = 'PDF is empty'
    text_chunks = read_and_split_txt(txt_contents, max_tokens)

    results = {}
    prompt_tokens = 0
    completion_tokens = 0
    for chunk in text_chunks:
        payload = {
            "model": 'deepseek-v3',
            "max_tokens": max_tokens,
            "temperature": temperature,
            "messages": [
                {
                    "role": "user",
                            "content":
                            # f"""
                            # Given a chemical reaction formula: {chemical_reaction}, please find the chemical name corresponding to each chemical formula from the following PDF text.
                            # The PDF text is: {chunk}.
                            # Search requirements are as follows:
                            # First, extract all chemical molecular formulas from the chemical reaction equation, then search only for the chemical names corresponding to the molecular formulas that appear in this reaction equation, Output as many molecular formulas as there are molecules, do not include any extra output.
                            # The output requirements are as follows:
                            # 1. Please output the results strictly according to "Chemical molecular Formula: Chemical Name".
                            # 2. Convert the chemical formula to a standard chemical reaction formula format. Make sure the output does not use Unicode subscript characters, LaTeX syntax, or any special symbols.
                            # """
                            f"""
                            Given a chemical reaction formula: {chemical_reaction}, please search the PDF text for the chemical names corresponding to each reactant chemical formula and product chemical formula in this chemical reaction (e.g., the chemical name corresponding to 'H2O' is 'water'). If the PDF text is empty, simply output according to the requirements without providing any explanatory text.
                            The PDF text is: {chunk}.

                            Search requirements are as follows:
                            1. First, extract all chemical formulas in the chemical reaction formula, and only search for the chemical names corresponding to the chemical formulas that appear in this chemical reaction.

                            Output requirements are as follows:
                            1. The output must be in English, and the results must strictly follow the format "chemical formula: chemical name" or "chemical formula: None", without any additional explanatory text.

                            2. Convert the chemical formulas into the standard chemical reaction format. Ensure the output does not use Unicode subscript characters, LaTeX syntax, or any special symbols.
                            """
                }
            ]
        }
        response = make_post_request(url, headers, payload)

        if response != None and response.status_code == 200:
            data = json.loads(response.text)
            if 'choices' in data:
                result_text = data['choices'][0]['message']['content']
                prompt_tokens += data['usage']['prompt_tokens']
                completion_tokens += data['usage']['completion_tokens']

                lines = result_text.split('\n')

                for line in lines:
                    if has_chinese(line) or ':' not in line:
                        continue

                    try:
                        molecule, name = line.split(': ', 1)
                    except:
                        continue
                    if '$' in name or len(name) <= 3:
                        name = 'None'
                    if molecule in results:
                        if name != 'None':
                            results[molecule] = name
                    else:
                        results[molecule] = name
        else:
            print("Request failed!")
            continue

    return results, prompt_tokens, completion_tokens
