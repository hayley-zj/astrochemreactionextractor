import argparse
import os
import json
import math
import numpy as np

def get_file_path_list(parent_dir: str, file_list: list) -> list:
    if len(file_list) == 0:
        return file_list

    file_path_list = []
    for file in file_list:
        if os.path.exists(file):
            file_path_list.append(file)
            continue

        if file[0] == '/':
            file_path = file
        else:
            file_path = os.path.join(parent_dir, file)
        if os.path.exists(file_path):
            file_path_list.append(file_path)

    return file_path_list


class ParamController:
    def __init__(self, desc='Process astronomy algorithm.', usage=''):
        self.args = None
        self.parser = argparse.ArgumentParser(description=desc, usage=usage)
        self.parser.add_argument('--input_path', type=str, default='./papers_for_extraction', help='papers for extraction')
        self.parser.add_argument('--output_path', type=str, default='./output', help='output path')
        self.parser.add_argument('--query_path', type=str, default='./cache', help='query path')

    def get_arg_parser(self):
        return self.parser

    def get_args(self):
        self.args = self.parser.parse_args()
            
        if not os.path.exists(self.args.output_path):
            os.makedirs(self.args.output_path, 0o755)

        if not os.path.exists(self.args.query_path):
            os.makedirs(self.args.query_path, 0o755)

        return self.args

    def _get_file_path_list_from_input_params(self) -> list:
        file_list = []
        if self.args.input_files is not None and len(self.args.input_files) > 0:
            origin_file_list = self.args.input_files
            for file_name in origin_file_list:
                if os.path.exists(os.path.join(self.args.input_path, file_name)):
                    file_list.append(file_name)
        return file_list

    def query_file_path_list(self) -> list:
        file_list = self._get_file_path_list_from_input_params()
        return get_file_path_list(self.args.input_path, file_list)

    def get_output_file_name_suffix(self) -> str:
        return self.args.output_name_suffix

    def save_output_file_list(self, file_name: str):
        output_list_file = os.path.join(self.args.log_path, 'output_files.list')
        with open(output_list_file, 'a') as f:
            f.write(file_name + '\n')

    def save_result_list(self, results: list):
        # 清理结果中的 NaN 和 Inf 值
        cleaned_results = self._cleanse_results(results)
        
        result_file = os.path.join(self.args.log_path, 'result_list.json')
        with open(result_file, 'w') as f:
            json.dump(cleaned_results, f, cls=NumpyEncoder, indent=4)

    def _cleanse_results(self, results: list) -> list:
        cleaned_results = []
        for item in results:
            if isinstance(item, dict):
                cleaned_item = {}
                for k, v in item.items():
                    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
                        cleaned_item[k] = None
                    else:
                        cleaned_item[k] = v
                cleaned_results.append(cleaned_item)
            else:
                cleaned_results.append(item)
        return cleaned_results


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        return super(NumpyEncoder, self).default(obj)