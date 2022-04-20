import requests

def parse(build):
    # set containing all possible files
    fileList={}
    url = "hgdownload.soe.ucsc.edu/goldenPath/"+ build+"/liftOver/"
    # parsing source code
    for line in str(requests.get("https://"+url).content).split("\\n"):
        line=line.strip()
        # detecting possible file link
        if "<a href=\"hg" in line:
            fileName=""
            start = False
            # extracting
            for character in line:
                if start == True and character == "<": # end file name
                    break
                elif start==True: # reading file name
                    fileName += character
                elif character ==">": # start parsing file name
                    start = True
            # adding file name with url download link
            fileNameconv = fileName.removesuffix('.over.chain.gz')
            fileNameconv = fileNameconv[6:]
            fileList[fileNameconv]=url+fileName
    return fileList