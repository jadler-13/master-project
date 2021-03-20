#changes space group number in SIMULATE.INP in "ref" and "new" directories;

with open("spacelog", "r") as fin:
    for line in fin:
        lc = line.split()
        if not lc:
            continue
        elif lc[0] == "*" and lc[1] == "Space" and lc[2] == "group":
            sg_number = lc[-1][:-1]
            break
        else:
            continue

with open("SIMULATE.INP", "r") as fin:
    with open("temp_SIMULATE.INP", "w") as fout:
        for line in fin:
            lc = line.split()
            if not lc:
                fout.write(line)
            elif lc[0] == "SPACE_GROUP_NUMBER=":
                fout.write(line.replace(lc[1], sg_number))
            else:
                fout.write(line)

with open("../new/SIMULATE.INP", "r") as fin:
    with open("../new/temp_SIMULATE.INP", "w") as fout:
        for line in fin:
            lc = line.split()
            if not lc:
                fout.write(line)
            elif lc[0] == "SPACE_GROUP_NUMBER=":
                fout.write(line.replace(lc[1], sg_number))
            else:
                fout.write(line)

