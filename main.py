# to run: python main.py build
# example: python main.py hg19

from email.headerregistry import ContentTransferEncodingHeader
import urlParser

import sys
import os



def main(argv):
    fileList= urlParser.parse(sys.argv[1])
    running = True
    print("-------------------------------------------------------------------------------")
    print("Please enter the chain file you would like to download from the available files")
    print("-------------------------------------------------------------------------------")
    for i in fileList:
        print(i)
    while(running):
        file = input("Type conversion (i.e. AilMel1 ). Type \"Quit\" at anytime to cancel: " ).strip()
        if file == "Quit":
            running = False
        elif file not in fileList:
            print("File not found. Try again.")
        else:
            userInput = input("Type \"Yes\" to confirm download or press any key to select another file: ")
            if userInput== "Yes":
                print("Downloading "+fileList[file])
                command = "wget --timestamping \'ftp://"+fileList[file] + "\' -O library/chainfiles/"+sys.argv[1]+"To"+file+".over.chain.gz"
                print(command)
                os.system(command)
                userInput = input("Type \"Continue\" to download another file or press any key to continue to triple lift over: ")
                if userInput == "Continue":
                    continue
                else:
                    running = False
            elif userInput == "Quit":
                running = False
            else:
                continue

if __name__ == "__main__":
   main(sys.argv[1:])