import sys

input_file = sys.argv[1]

with open(input_file) as fh:
    header = next(fh)  # 跳过表头

    start_state = 0
    state_con = 0
    prev_line = None

    for line in fh:
        line = line.rstrip("\n")
        temp = line.split(",")

        state = int(temp[4])

        if state == start_state:
            state_con += 1
        else:
            if state_con >= 5 and start_state != 0 and prev_line is not None:
                # 使用上一行作为区间终点
                end = int(prev_line[2])
                start = end - state_con
                chrom = prev_line[1]
                print(f"{chrom}\t{start}\t{end}\t{start_state}")

            start_state = state
            state_con = 1

        prev_line = temp

    # ===== 处理文件结尾（Perl 原代码这里是漏掉的）=====
    if state_con >= 5 and start_state != 0 and prev_line is not None:
        end = int(prev_line[2])
        start = end - state_con
        chrom = prev_line[1]
        print(f"{chrom}\t{start}\t{end}\t{start_state}")
