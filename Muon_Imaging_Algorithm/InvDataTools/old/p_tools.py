# author:高金磊
# datetime:2021/11/29 12:38
if __name__ == '__main__':
    file=open("E:\vscode\Muon_Imaging_Algorithm\dataTools\data\p",'r')
    res=[]
    while 1:
        line=file.readline().replace(' ','').replace("\n",'')
        if line:
            try:
                res.append(float(line))
            except Exception:
                pass
        else:
            break
    file.close()
    file=open("E:\vscode\Muon_Imaging_Algorithm\dataTools\data\p",'w')
    middle=1
    for re in res:
        file.write(str(middle/10000))
        file.write('\n')
        middle+=1
    file.close()
