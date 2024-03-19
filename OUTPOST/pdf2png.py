# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 15:36:07 2024

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

================================ description ==================================
本来是使用pdf2image处理的，但是linux下要安装poppler，且要sudu权限（conda装不了，因为OUTPOST
环境有很多是pip安装的，冲突），所以现在用fitz处理。原来的pdf2image脚本注释掉了，在windows环境下
也可以用。
=================================== input =====================================

=================================== output ====================================

================================= parameters ==================================

=================================== example ===================================

=================================== warning ===================================
"""
import fitz  # 又名PyMuPDF
import sys
import os

def convert_pdf_last_page_to_png(pdf_dp, png_dp, initial_dpi=600):
    # 打开PDF文档
    doc = fitz.open(pdf_dp)
    # 获取最后一页
    last_page = doc[-1]
    
    # 定义渲染图像并检查大小的函数
    def render_and_check_size(page, dpi, output_path):
        zoom_factor = dpi / 72
        matrix = fitz.Matrix(zoom_factor, zoom_factor)
        pix = page.get_pixmap(matrix=matrix)
        pix.save(output_path)
        return os.path.getsize(output_path)
    
    # 首先使用初始DPI渲染图像
    file_size = render_and_check_size(last_page, initial_dpi, png_dp)
    
    # 如果文件大小超过2MB（2 * 1024 * 1024字节），使用72 DPI重新渲染
    if file_size > 2 * 1024 * 1024:
        render_and_check_size(last_page, 72, png_dp)
    
    # 关闭文档
    doc.close()

if __name__ == '__main__':
    # 接收命令行参数指定的目录
    dir_ = sys.argv[1]
    # dir_ = r"D:\Projects\GEMINI\outpost2\figs"
    pdfs = os.listdir(dir_)
    for pdf in pdfs:
        if pdf.endswith('.pdf'):
            pdf_dp = os.path.join(dir_, pdf)
            png_dp = pdf_dp.replace('.pdf', '.png')
            convert_pdf_last_page_to_png(pdf_dp, png_dp)




# import pdf2image
# import sys
# import os

# def convert_pdf_last_page_to_png(pdf_dp, png_dp):
#     # 将PDF文档的最后一页转换成图像
#     images = pdf2image.convert_from_path(pdf_dp)
#     last_page = images[-1]  # 获取最后一页
#     # 保存最后一页为PNG文件
#     last_page.save(png_dp, 'PNG')
    

# if __name__ == '__main__':
    
#     # dir_ = r"D:\Projects\GEMINI\results2_cp\biomarkers_analysis\ANCOM_identification"
#     dir_ = sys.argv[1]
#     pdfs = os.listdir(dir_)
#     for pdf in pdfs:
#         if pdf.endswith('.pdf'):
#             pdf_dp = os.path.join(dir_, pdf)
#             png_dp = pdf_dp.replace('.pdf', '.png')
#             convert_pdf_last_page_to_png(pdf_dp, png_dp)
    
#     # pdf_dp = sys.argv[1]
#     # png_dp = sys.argv[2]
    
#     # convert_pdf_last_page_to_png(pdf_dp, png_dp)
    
    
    


