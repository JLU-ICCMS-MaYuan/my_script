import re

def read_train_cfg(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    return content

def write_train_cfg(file_path, content):
    with open(file_path, 'w') as file:
        file.write(content)

def filter_cfg_data(content, min_distance_threshold=0.5):
    # 匹配BEGIN_CFG到END_CFG之间的所有内容
    pattern = re.compile(r'BEGIN_CFG.*?END_CFG', re.DOTALL)
    cfg_blocks = pattern.findall(content)

    filtered_blocks = []
    for block in cfg_blocks:
        # 匹配mindist的值
        mindist_match = re.search(r'Feature\s+mindist\s+([\d.]+)', block)
        if mindist_match:
            mindist_value = float(mindist_match.group(1))
            # 如果mindist值大于等于阈值，则保留该块
            if mindist_value >= min_distance_threshold:
                filtered_blocks.append(block)
    
    # 将过滤后的块重新组合成一个字符串
    return '\n'.join(filtered_blocks)

def main():
    input_file = 'train.cfg'
    output_file = 'filtered_train.cfg'
    
    # 读取train.cfg文件内容
    content = read_train_cfg(input_file)
    
    # 过滤数据块
    filtered_content = filter_cfg_data(content)
    
    # 将过滤后的内容写入新的文件
    write_train_cfg(output_file, filtered_content)
    print(f"Filtered content written to {output_file}")

if __name__ == "__main__":
    main()
