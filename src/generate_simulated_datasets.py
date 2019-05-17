import json
import simulate

def main(configs):
    raise NotImplementedError
    
if __name__ == '__main__':
    input_json = '../data/external/experiments.json'
    with open(input_json) as f:
        configs = json.load(f)