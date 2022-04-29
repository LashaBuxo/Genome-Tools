import sys

file = open('test.fa', 'r')
Lines = file.readlines()

sequence = ""
for line in Lines:
    sequence += line.replace('\n','')


max_repeat = ""
for i in range(0, len(sequence)):
    for repeat_length in range(1, 6):
        repeat_short_seq = sequence[i:(i + repeat_length)]
        index = i
        build_seq = repeat_short_seq
        while index < len(sequence):
            index += repeat_length
            cur_short_seq = sequence[index:(index + repeat_length)]
            if cur_short_seq != repeat_short_seq: break
            build_seq += repeat_short_seq
        if len(build_seq) > len(max_repeat):
            max_repeat = build_seq

print(max_repeat)
