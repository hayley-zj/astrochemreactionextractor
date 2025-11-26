
import time
from astrochemreactionextractor import logger
from astrochemreactionextractor.utils.params import ParamController
from astrochemreactionextractor.utils.utils import wrap_with_exception_handler
from astrochemreactionextractor.reaction_extraction_pipe import ReactionPipeline

@wrap_with_exception_handler
def init():
    pi = ParamController()
    args = pi.get_args()

    return args.input_path, args.output_path, pi

if __name__ == '__main__':
    # This chemical reaction extraction algorithm example can be run on server XXXX.
    # python main.py --input_path *** --output_path ***

    input_path, output_path, pi = init()

    # Extraction started
    pipline  = ReactionPipeline(input_path, output_path)
    pipline.extract()