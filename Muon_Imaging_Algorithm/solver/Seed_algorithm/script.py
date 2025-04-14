# author:高金磊
# datetime:2022/9/15 10:00
from DataTools.Jxyz_Tools import getxyz_from_shape


def create_refs_bounds():
    path=r"M:\pycharm\Inversion\data\seed_study"
    refs_file=open(path+r"\refs",'w')
    bounds_file=open(path+r"\bounds","w")
    # for i in range(82*82*82):
    #     x,y,z=getxyz_from_shape((82,82,82),1+i)
    #     if y == 1 :
    #         refs_file.write('0.0')
    #         refs_file.write("\n")
    #         bounds_file.write("-0.1 0.1")
    #         bounds_file.write("\n")
    #         continue
    #     refs_file.write('2.65')
    #     refs_file.write("\n")
    #     bounds_file.write("0 2.7")
    #     bounds_file.write("\n")

    for z in range(64):
        for y in range(64):
            for x in range(64):
                if x==0 or y==0 or z==0 or x==63 or y==63 or z==63:
                    refs_file.write('0.0')
                    refs_file.write("\n")
                    bounds_file.write("-0.1 0.1")
                    bounds_file.write("\n")
                    continue
                if x==1 or y==1 or z==1 or x==62 or y==62 or z==62:
                    refs_file.write('0.0')
                    refs_file.write("\n")
                    bounds_file.write("-0.1 0.1")
                    bounds_file.write("\n")
                    continue
                if x==2 or y==2 or z==2 or x==61 or y==61 or z==61:
                    refs_file.write('2.65')
                    refs_file.write("\n")
                    bounds_file.write("2.64 2.67")
                    bounds_file.write("\n")
                    continue
                refs_file.write('2.65')
                refs_file.write("\n")
                bounds_file.write("0 2.7")
                bounds_file.write("\n")
    refs_file.close()
    bounds_file.close()
if __name__ == '__main__':
    create_refs_bounds()