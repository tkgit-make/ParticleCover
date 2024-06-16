def main():
    filename = "cmake-build-debug/differenceFinal.txt"
    with open(filename) as file:
        entriesList = [] 
        currentEntry = [] 
        for line in file:
            brokenLine = line.rstrip().split()

            if line.find("c") != -1 or line.find("a") != -1 or line.find("d") != -1:
                entriesList.append(currentEntry)
                currentEntry = []
                currentEntry.append(brokenLine)
            else:
                currentEntry.append(brokenLine)

    for entry in entriesList:
        if ['---'] not in entry:
            print(entry)
            continue

        inputNums = []
        outputNums = [] 
        for j in range(1, len(entry)):
            if entry[j][0] == "<":
                for i in range(1, len(entry[j])):
                    if entry[j][i][0] == "[":
                        inputNums.append(float(entry[j][i][1:-1]))
                    elif entry[j][i][-1] == "]":
                        inputNums.append(float(entry[j][i][0:-1]))
                    else: 
                        if entry[j][i] == 'Superpoint':
                            print(entry)
                            print("---")
                        else:
                            inputNums.append(float(entry[j][i]))

            if entry[j][0] == ">":
                for i in range(1, len(entry[j])):
                    if entry[j][i][0] == "[":
                        outputNums.append(float(entry[j][i][1:-1]))
                    elif entry[j][i][-1] == "]":
                        outputNums.append(float(entry[j][i][0:-1]))
                    else: 
                        if entry[j][i] == 'Superpoint':
                            print(entry)
                            print("---")
                        else:
                            outputNums.append(float(entry[j][i]))

        diffThreshold = 0.001
        multiple = 10000
        if len(inputNums) != len(outputNums):
            print(entry)
            print("---")
            continue

        for i in (range(len(inputNums))):
            if abs(inputNums[i]) > 1000: 
                if abs(inputNums[i] - outputNums[i]) > diffThreshold * multiple:
                    #print("input", inputNums)
                    #print("output", outputNums)
                    print(entry)
                    print("---")
                    break 
            else: 
                if abs(inputNums[i] - outputNums[i]) > diffThreshold:
                    #print("input", inputNums)
                    #print("output", outputNums)
                    print(entry)
                    print("---")
                    break 

        
        

        
        
        

if __name__ == "__main__":
    main()