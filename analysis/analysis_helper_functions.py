


# iterate a file linewise
# check line for keyword
# if keyword is found return value after seperator
def fileValueFromKeyword(filepath, keyword, seperator='='):
    found_counter = 0
    assert(filepath)
    value = None
    try:
        with open(filepath) as FILE:
            for line in FILE:
                if keyword in line:
                    found_counter += 1
                    value = line.split(seperator,1)[1].rstrip()
                    break
    except:
        pass
    if found_counter == 0:
        raise Exception("unable to find keyword "+keyword+" in file "+filepath+"\n")
    # if found_counter >= 2:
        # raise Exception("multiple keywords "+keyword+" in file "+filepath+"\n")
    return value


