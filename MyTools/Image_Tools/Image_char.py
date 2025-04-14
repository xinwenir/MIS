# author:高金磊
# datetime:2022/8/1 14:48

from PIL import Image
def _getchar(r, g, b):
    #####彩图
    # txt = list("$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/\|()1{}[]?-_+~<>i!lI;:,\"^`'.")
    ####黑白图
    txt = list("* ")
    gray = 0.2126 * r + 0.7152 * g + 0.0722 * b
    bili = 256 / len(txt)
    index = int(gray / bili)
    return txt[index]
def get_txt_image(photo="M:\pycharm\pythonProject\Tools\MyTools\Char_image\签名.png"):
    image = Image.open(photo)  # 导入一张图片
    width = 107
    height = 38
    image = image.resize((width, height), Image.NEAREST)
    txt = ""
    for y in range(height):
        for x in range(width):
            pix = image.getpixel((x, y))
            txt += _getchar(r=pix[0],g=pix[1],b=pix[2])
        txt += "\n"
    return txt
# print(get_txt_image())