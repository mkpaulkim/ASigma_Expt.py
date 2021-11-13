from tkinter import filedialog
from tkinter import messagebox
import cv2
import numpy as np
import pycubelib.data_functions as df
import pycubelib.general_functions as gf
from os import path

u16limit = (0, 2 ** 16 - 1)


def write_txt(text, txt_path=""):
    if not txt_path:
        txt_path = filedialog.asksaveasfilename(title='.TXT file name to write',
                                                filetypes=[('TXT files', '*.txt *.TXT')],
                                                defaultextension='.txt')
    if path.exists(txt_path):
        ans = messagebox.askyesno('file exists', '.txt file exists - overwrite? ')
        if not ans:
            return ''

    with open(txt_path, 'w') as f:
        f.write(text)
    return txt_path


def read_txt(txt_path=''):
    if not txt_path:
        txt_path = filedialog.askopenfilename(title='open .TXT file',
                                              filetypes=[('TXT files', '*.txt *.TXT')])
    with open(txt_path) as f:
        text = f.read()

    return text, txt_path


def write_png(aa, png_path='', alimit=()):
    """
    write to 16 bit gray png file
    """

    if not png_path:
        png_path = filedialog.asksaveasfilename(title='.PNG file name to write',
                                                filetypes=[('PNG files', '*.png *.PNG')],
                                                defaultextension='.png')
    if alimit:
        aa = df.renormalize(aa, alimit, u16limit, np.uint16)
    else:
        alimit = (np.min(aa), np.max(aa))
        aa = aa.astype(np.uint16)

    cv2.imwrite(png_path, aa)

    return alimit, png_path


def read_png(png_path='', alimit=()):
    """
    read from 16 bit gray png file
    """

    if not png_path:
        png_path = filedialog.askopenfilename(title='.PNG file to read',
                                              filetypes=[('PNG files', '*.png *.PNG')])
    aa = cv2.imread(png_path, cv2.IMREAD_UNCHANGED)
    if alimit:
        aa = df.renormalize(aa, u16limit, alimit, type(alimit[0]))
    else:
        alimit = (np.min(aa), np.max(aa))
    nyx = np.shape(aa)

    file_info = {'nys':nyx, 'alimit':alimit, 'png_path':png_path}
    return aa, file_info


if __name__ == '__main__':
    aa, _, _, _ = read_png()
    gf.what_is('aa', aa)

    _, path_b = write_png(aa)

    bb, _, _, _ = read_png(path_b)

    gf.what_is('bb', bb)

    text, _ = read_txt()
    print(text)
    write_txt(text)



