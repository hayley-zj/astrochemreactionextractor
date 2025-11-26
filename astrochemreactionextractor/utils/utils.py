import traceback
from functools import wraps
import re
import time
import fitz  # PyMuPDF
import requests
import pandas as pd

# from utils.logger import logger


def wrap_with_exception_handler(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            # logger.error("{} exception: {}".format(func.__name__, traceback.format_exc()))
            exit(-1)
    return wrapper


# Define the chemical element symbols in the periodic table
periodic_table_elements = {
    "c", "s", "E", "l", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
    "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs"
}

witelist_elements = {'photon', 'Photon', 'e', 'CR', 'CRPHOT', 'CRP', 'hu', 'hnu', 'h'}

def need_verify_validity_element(element):
    """
    判断元素符号的长度是否小于2且不含数字

    :param element: str, 元素符号
    :return: bool, 如果元素符号长度小于2且不含数字返回True，否则返回False
    """
    # 判断长度小于2
    if len(element) >= 2:
        return False

    # 判断是否包含数字
    if any(char.isdigit() for char in element):
        return False

    return True

def is_valid_formula(formula_list):
    for formula in formula_list:
        try:
            # formula = formula.encode().decode('unicode_escape')
            formula = formula.encode('latin1').decode('utf-8')
        except (UnicodeEncodeError, UnicodeDecodeError) as e:
            formula = formula.encode('utf-8').decode('utf-8')
        # 使用正则表达式过滤非英文字母字符
        formula = re.sub(r'\(.*?\)', '', formula) 
        formula = re.sub(r'[^a-zA-Z]', '', formula)
        if formula in witelist_elements:
            continue
        length = len(formula)
        if length <= 0:
            return False
        
        i = 0
        while i < length:
            if i == length - 1:  # 最后一个字母
                if formula[i] not in periodic_table_elements:
                    return False
                i += 1
            else:
                if formula[i:i+2] in periodic_table_elements:
                    i += 2
                elif formula[i] in periodic_table_elements:
                    i += 1
                else:
                    return False
    return True

def check_elements_in_periodic_table(formula_list):
    """
    判断formula_list中的每一个元素符号是否都在化学周期表中

    :param formula_list: list, 元素符号的列表
    :return: bool, 如果所有元素符号都在周期表中返回True，如果有一个不在则返回None
    """
    for element in formula_list:
        if need_verify_validity_element(element) and element not in periodic_table_elements:
            return False
    return True

def get_smiles_from_simple_species(formula):
        return 'None'

def clean_text(text):
    # Remove line breaks and join split words caused by hyphens at line endings
    text = re.sub(r'-\n', '', text)
    # text = re.sub(r'\n', ' ', text)
    return text

def extract_doi_based_on_condition(input_string):
    if re.match(r'^\d', input_string):  # 判断字符串是否以数字开头
        # 提取前两个 '.' 的内容
        parts = input_string.split('.', 2)
        if len(parts) > 2:
            return f"{parts[0]}.{parts[1]}"
        else:
            return input_string  # 如果只有一个 '.'，返回整个字符串
    else:
        # 提取第一个 '.' 之前的内容
        parts = input_string.split('.', 1)
        return parts[0]

def get_smiles_from_pubchem(chemical_name, max_retries=4, retry_delay=2):
    retries = 0
    while retries < max_retries:
        # 构建查询URL
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/property/CanonicalSMILES/JSON"
        
        # 发送请求
        try:
            response = requests.get(url)
            
            # 检查响应状态
            if response.status_code == 200:
                # 解析JSON响应
                data = response.json()
                try:
                    # 提取SMILES字符串
                    smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
                    return smiles
                except (KeyError, IndexError):
                    return "Error: Unable to find SMILES for the given chemical name."
            else:
                return f"Error: API request failed with status code {response.status_code}."
        except (requests.exceptions.ConnectionError, requests.exceptions.ChunkedEncodingError) as e:
            print(f"ConnectionError occurred: {e}")
            print(f"Retrying in {retry_delay} seconds...")
            retries += 1
            time.sleep(retry_delay)
    return 'None'

def extract_reactions_from_tex(tex_file_path):
    # 打开并读取.tex文件
    with open(tex_file_path, 'r', encoding='utf-8') as file:
        tex_content = file.read()
    
    # 定义正则表达式
    reaction_pattern = re.compile(
        r'([A-Za-z$_0-9\\\^\{\}\(\)\-\+\s]+\s*\\rightarrow\s*[A-Za-z$_0-9\\\^\{\}\(\)\-\+\s]+)'
        # r'|([A-Za-z$_0-9\(\)\-\+\s]+\s*\\rightarrow\s*[A-Za-z$_0-9\(\)\-\+\s]+)'
    )
    equations = reaction_pattern.findall(tex_content)
    
    # 将匹配的结果格式化为化学反应式
    formatted_equations = [f"{eq[0].strip()} -> {eq[1].strip()}" for eq in equations]

    retrieval_context_txt = '\n'.join([item for item in formatted_equations])
    
    return retrieval_context_txt

def extract_title_from_tex(tex_content):
    reaction_pattern = re.compile(
        r'^.*\\title.*$', re.MULTILINE
    )
    # Find all matches in the .tex content
    titles = reaction_pattern.findall(tex_content)
    
    if titles:
        return titles[0]
    else:
        return None

def extract_reactions_txt_from_tex(tex_content):
    reaction_pattern = re.compile(
        r'(\\ce\{[^}]+\}(?:\s*\+\s*\\ce\{[^}]+\})*\s*&\s*\\hspace\{[^}]+\}\$\s*\\longrightarrow\$\s*\\ce\{[^}]+\}(?:\s*\+\s*\\ce\{[^}]+\})*)'
        r'|(&\s*\\hspace\{[^}]+\}\$\s*\\longrightarrow\$\s*\\ce\{[^}]+\}(?:\s*\+\s*\\ce\{[^}]+\})*)'
        r'|([A-Za-z$&~*_0-9\\.\^\{\}\(\)\-\+\s]{1,}\s*(\\rightarrow|\\longrightarrow|\\to|->|\$\s*\\to\$)\s*[A-Za-z$&~*_0-9\\.\^\{\}\(\)\-\+\s]+)'
        r'|\&\s*~+\s*(.*?)\s*\&'
    )

    # Find all matches in the .tex content
    reactions = reaction_pattern.findall(tex_content)
    formatted_equations = [next((part for part in reaction if part), '') for reaction in reactions]

    retrieval_context_txt = '\n'.join([item for item in formatted_equations])
    return retrieval_context_txt

def read_and_split_txt(txt_content, max_tokens):   
    # Split content by lines
    lines = txt_content.split('\n')
    
    chunks = []
    current_chunk = ""
    
    for line in lines:
        if len(current_chunk) + len(line) + 1 <= max_tokens:
            # Add the line to the current chunk
            if current_chunk:
                current_chunk += '\n'
            current_chunk += line
        else:
            # Add the current chunk to chunks and start a new chunk
            chunks.append(current_chunk)
            current_chunk = line
    
    # Add the last chunk if it's not empty
    if current_chunk:
        chunks.append(current_chunk)
    
    return chunks


def extract_formula_from_recation(reaction):
    # 去除空格并分离反应物和生成物
    # reaction = reaction.replace(' ', '')

    # 使用正则表达式匹配化学分子式，包括下标和上标
    pattern = re.compile(
        r'\\ce{([^}]+)}' 
        r'|(\$?[a-z]?\$?-?[A-Z][a-z]?(?:\$_\d+\$)?(?:[A-Z][a-z]?(?:\$_\d+\$)?)*)'
        r'|([A-Za-z][A-Za-z0-9_{}^+-]*)'
    )

    # 提取反应物和生成物
    formula_list = pattern.findall(reaction)
    formatted_formula_list = [next((part for part in formula if part), '') for formula in formula_list]

    return formatted_formula_list


def find_formulas_contexts(tex_content, chemical_list, context_size=8):
    contexts = []
    # Split the text into words
    words = tex_content.split()
    for chemical in chemical_list:
        # 使用 re.escape 转义 chemical 以处理特殊字符
        escaped_chemical = re.escape(chemical)
        matches = []
        try: 
            # pattern = re.compile(r'\([^)]*' + escaped_chemical + r'[^)]*\)')
            pattern = re.compile(escaped_chemical)
            # pattern = re.compile(r'\b' + re.escape(escaped_chemical) + r'\b')
        except re.error as e:
            # 如果正则表达式有错误，输出错误信息但不停止程序
            print(f"正则表达式错误: {e}")
            return ''
        # 查找所有匹配的括号内容
        matche_str = [match.group() for match in pattern.finditer(tex_content)]
        if len(matche_str) != 0:
            # matches = [m.start() for m in re.finditer(matche_str[0], text)]
            matches = [m.start() for m in re.finditer(re.escape(matche_str[0]), tex_content)]
        else:
            matches = []

        for match in matches:
            # Find the word index where the match occurs
            word_index = len(re.findall(r'\S+', tex_content[:match]))
            # Extract context around the match
            start = max(word_index - context_size, 0)
            end = min(word_index + 0, len(words))
            context = ' '.join(words[start:end])
            contexts.append(context)
    retrieval_contexts = '\n'.join([item for item in contexts])
    return retrieval_contexts

def extract_contents_reactions_from_pdf(pdf_document):
    document = fitz.open(pdf_document)

    pdf_reactions = ''
    pdf_contents = ''

    for page_num in range(document.page_count):
        page = document.load_page(page_num)
        text = page.get_text()
        pdf_contents += text
        print('page_num',page_num)

        ####定义一个新的正则表达式模式
        pattern = re.compile(
            r'([\w\s\−\,\·\⋅\.\*\∗()\+\-]+\s*→\s*[\w\s\−\,\·\⋅\.\*()\+\-]+(?:\s*[+þ]\s*[\w\s\−\,\·\⋅\.\*\∗()\+\-]+)*)'
            r'|(→\s*[\w\s\−\,\·\⋅\.\*\∗()\+\-]+(?:\s*\+\s*[\w\s\−\,\·\⋅\.\*\∗()\+\-]+)*)'
        )
        # 查找所有匹配项
        matches = pattern.findall(text)
        if len(matches) != 0:
            formatted_equations = [next((part for part in reaction if part), '') for reaction in matches]

            pdf_reactions += '\n' + '\n'.join([item for item in formatted_equations])
    return pdf_contents, pdf_reactions


def extract_match_coordinates_from_blocks(pdf_path):
    # 定义正则表达式模式
    pattern = re.compile(
        r'([\w\s\−\,\·\⋅\.\*\∗()\+\-]+\s*→\s*[\w\s\−\,\·\⋅\.\*()\+\-]+(?:\s*[+þ]\s*[\w\s\−\,\·\⋅\.\*\∗()\+\-]+)*)'
        r'|(→\s*[\w\s\−\,\·\⋅\.\*\∗()\+\-]+(?:\s*\+\s*[\w\s\−\,\·\⋅\.\*\∗()\+\-]+)*)'
    )

    # 打开PDF文件
    doc = fitz.open(pdf_path)
    pdf_reactions = ''
    pdf_contents = ''
    
    data = []
    # 遍历每一页
    for page_num in range(len(doc)):
        page = doc.load_page(page_num)
        pdf_contents += page.get_text()
        
        # 提取文本块
        blocks = page.get_text("blocks")
        
        # 遍历文本块
        for block in blocks:
            text = block[4].strip()  # 文本内容
            bbox = fitz.Rect(block[:4])  # 文本块的边界框
            
            # 使用正则表达式匹配文本
            matches = pattern.findall(text)
            if matches:
                formatted_equations = [next((part for part in reaction if part), '') for reaction in matches]
                pdf_reactions += '\n' + '\n'.join([item for item in formatted_equations])
                data.append([formatted_equations, bbox])
                # print(f"===== Page {page_num + 1} =====")
                # print(f"Match: {matches}")
                # print(f"Coordinates (x0, y0, x1, y1): {bbox}")
                # print("------")
    df = pd.DataFrame(data, columns=['reactions_txt', 'coordinates'])
    return df