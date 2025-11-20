import pyfiglet

def qe_logo(text, font="standard"):
    """
    使用 pyfiglet 生成大号的 ASCII 文字
    :param text: 要转换的文本
    :param font: 使用的字体（默认为 'standard'）
    """
    figlet = pyfiglet.Figlet(font=font)
    big_text = figlet.renderText(text)
    print(big_text)

if __name__ == "__main__":
    text = "SCRIPTS 4 QE"
    qe_logo(text)