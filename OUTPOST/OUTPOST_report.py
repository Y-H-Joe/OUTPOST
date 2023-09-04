from datetime import datetime
from pdf2image import convert_from_path
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
# from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Table, TableStyle
# import pandas as pd
import csv
import re
import sys

def update_csv_content(input_filename, output_filename, replace_dict):
    rows = []

    # 读取 CSV 文件
    with open(input_filename, 'r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        for row in reader:
            updated_row = []
            for cell in row:
                new_cell = cell
                for old_text, new_text in replace_dict.items():
                    new_cell = new_cell.replace(old_text, new_text)
                updated_row.append(new_cell)
            rows.append(updated_row)

    # 将更新后的数据写到新的 CSV 文件
    with open(output_filename, 'w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerows(rows)

# Functions
def _set_font(used_canvas, font_name, size):
    used_canvas.setFont(font_name, size)

def insert_header_footer(used_canvas, header_text, footer_text, is_cover=False):
    global current_page_number

    _set_font(used_canvas, "Helvetica", HEADER_FONT_SIZE)
    used_canvas.drawString(HEADER_X, HEADER_Y, header_text)

    _set_font(used_canvas, "Helvetica", FOOTER_FONT_SIZE)

    if is_cover:
        footer_text = f"{footer_text} {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    else:
        current_page_number += 1
        footer_text = f"Page: {current_page_number}"

    used_canvas.drawRightString(page_width - FOOTER_RIGHT_MARGIN, FOOTER_Y, footer_text)

def convert_pdf_page_to_image(pdf_path, page_number):
    images = convert_from_path(pdf_path, first_page=page_number, last_page=page_number)
    return images[0] if images else None

def insert_image_from_pdf(used_canvas, image_stream, x, y, width, height):
    """Insert an image onto the canvas at the specified position."""
    used_canvas.drawInlineImage(image_stream, x, y, width, height)

def insert_table_from_tsv(used_canvas, tsv_file, num_rows, num_cols, x, y, margin=20):
    # Adjust y position based on the top of the table
    y_top_of_table = y

    A4_WIDTH = A4[0]
    TARGET_TOTAL_WIDTH = A4_WIDTH - 2 * margin

    # Read data from tsv
    with open(tsv_file, "r", newline="") as file:
        raw_data = [line.strip().split("\t")[:num_cols] for line in file][:num_rows]

    # Define your styles
    style = getSampleStyleSheet()['BodyText']
    style.wordWrap = 'CJK'
    style.alignment = 1  # Center alignment

    data = []
    col_content_widths = [0] * num_cols

    for row_data in raw_data:
        processed_row = []
        for idx, cell_data in enumerate(row_data):
            cell_content = cell_data.replace("\n", "<br/>")
            para = Paragraph(cell_content, style)
            w, _ = para.wrap(TARGET_TOTAL_WIDTH, 10000)  # Temporarily set a large height to get natural dimensions

            col_content_widths[idx] = max(col_content_widths[idx], w)

            processed_row.append(para)
        data.append(processed_row)

    total_original_width = sum(col_content_widths)
    scaling_factor = TARGET_TOTAL_WIDTH / total_original_width

    col_widths = [w * scaling_factor for w in col_content_widths]

    max_heights = []
    for row_data in data:
        row_max_height = 0
        for idx, para in enumerate(row_data):
            _, h = para.wrap(col_widths[idx], 10000)  # 获取每个单元格的高度
            row_max_height = max(row_max_height, h + 15)  # +10 为了保留些许间距
        max_heights.append(row_max_height)

    # Calculate total table height
    total_table_height = sum(max_heights)

    # Adjust y position based on the total table height
    # y_adjusted = y - total_table_height

    # Define the table style
    table_style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.black),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 10),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 8),
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),
        ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 1), (-1, -1), 6),
        ('ALIGN', (0, 1), (-1, -1), 'LEFT'),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('TOPPADDING', (0, 0), (-1, -1), 4),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.black),
        ('BOX', (0, 0), (-1, -1), 0.5, colors.black),
    ])

    # Create table object
    table = Table(data, repeatRows=1, colWidths=col_widths, rowHeights=max_heights)
    table.setStyle(table_style)

    # Draw table with adjusted y position
    table.wrapOn(used_canvas, TARGET_TOTAL_WIDTH, total_table_height)
    table.drawOn(used_canvas, x, y_top_of_table - total_table_height)
    return total_table_height

def insert_text(used_canvas, x, y, text, font_size=12, max_width=None, bold=False):
    if bold:
        font = "Helvetica-Bold"
    else:
        font = "Helvetica"

    used_canvas.setFont(font, font_size)
    line_height = font_size + 4  # 增加行高以避免重叠

    if max_width is None:
        max_width = used_canvas._pagesize[0] - x

    # 先按照 ";" 拆分
    semi_colon_splits = text.split(';')

    lines = []

    for segment in semi_colon_splits:
        # 再按照空格、/或者\\拆分文本
        words = re.split(r'(\s|/|\\\\|_)', segment)
        words = [word for word in words if word]  # 去除空字符串

        current_line = ""

        for word in words:
            if word in ['/', '\\\\']:  # 如果是/或\\，则将它添加到最后一个单词上，而不是开启新的单词
                current_line += word
            else:
                test_line = current_line + "" + word if current_line else word
                line_width = used_canvas.stringWidth(test_line, font, font_size)

                if line_width <= max_width or word in ['/', '\\\\']:
                    current_line = test_line
                else:
                    lines.append(current_line)  # 将当前行添加到列表中
                    current_line = word

        if current_line:
            lines.append(current_line)

        # 如果不是最后一个segment, 加一个空行
        if segment != semi_colon_splits[-1]:
            lines.append('')

    for line in lines:
        y -= line_height  # 先减少行高，然后再绘制
        used_canvas.drawString(x, y, line)

    return y - 10  # 再次减少一些额外的空间以确保间隙

def set_cover(used_canvas, text="OUTPOST REPORT"):
    _set_font(used_canvas, "Helvetica", 20)
    x = page_width / 2
    y = page_height / 2
    used_canvas.drawCentredString(x, y, text)
    insert_header_footer(used_canvas, "OUTPOST", "Time:", True)  # is_cover set to True
    used_canvas.showPage()

def read_csv_content(file_path):
    with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        data = [row for row in reader]
    return data


def insert_and_update_y(used_canvas, y, path, description, tsv=False):
    # 打印路径并更新y
    if path:  # 检查路径是否为空
        y = insert_text(used_canvas, 20, y, f"Path: {path}", font_size=12, max_width=550)
        y -= 10  # update y for path text height

    # 打印解释内容并更新y
    y = insert_text(used_canvas, 20, y, description, font_size=12, max_width=550)  # 修改字体大小为10
    y -= 10  # Update y for description height

    if tsv and path:  # 路径不为空时才进行绘制
        y -= 20
        table_height = insert_table_from_tsv(used_canvas, path, 4, 7, 20, y)
        y -= table_height  # 使用表格的真实高度来更新y
    elif path:
        remaining_height = y - FOOTER_Y
        image_width = page_width - 2 * 20  # considering a 20 point margin on both sides
        image_height = (image_width * 3) / 4  # based on aspect ratio

        if image_height > remaining_height:
            image_height = remaining_height
            image_width = (image_height * 4) / 3

        image_x = (page_width - image_width) / 2  # 图片居中

        image_stream = convert_pdf_page_to_image(path, 1)
        insert_image_from_pdf(used_canvas, image_stream, image_x, y - image_height, image_width, image_height)
        y -= image_height  # Update y based on the image height
    return y

def create_pdf_content(content_csv, example_path_csv, used_canvas):
    content_data = read_csv_content(content_csv)
    example_data = read_csv_content(example_path_csv)

    for index, row in enumerate(content_data):
        insert_header_footer(used_canvas, "OUTPOST", "Page")
        y = page_height - 40  # slightly adjust y for better spacing
        y = insert_text(used_canvas, 20, y, row['title'], font_size=24, bold=True)

        # 检查title是否为"9. Notes"
        if row['title'] != "9. Notes":
            y = insert_text(used_canvas, 20, y, "Methods", font_size=14, max_width=550, bold=True)
            y = insert_text(used_canvas, 20, y, row['Methods'], font_size=12, max_width=550)
            y = insert_text(used_canvas, 20, y, "Results", font_size=14, max_width=550, bold=True)
            y = insert_text(used_canvas, 20, y, row['Results'], font_size=12, max_width=550)
        else:  # 当title为"9. Notes"时，直接输出Methods列中的内容，不输出Methods、Results标题
            y = insert_text(used_canvas, 20, y, row['Methods'], font_size=12, max_width=550)

        examples = [key for key in row.keys() if "Visualization example" in key]
        for idx, example in enumerate(examples):
            if row[example]:  # If there's content
                used_canvas.showPage()
                insert_header_footer(used_canvas, "OUTPOST", "Page")
                y = page_height - 40  # reset y to top

                y = insert_text(used_canvas, 20, y, example, font_size=14, max_width=550, bold=True)  # 修改字体大小为12
                path = example_data[index][example]
                description = row[example]
                if path:
                    y = insert_and_update_y(used_canvas, y, path, description, tsv=path.endswith('.tsv'))
                else:
                    y = insert_text(used_canvas, 20, y, description, font_size=12, max_width=550)

        used_canvas.showPage()


def generate_pdf(output_filename, explanation_upd, visualization_upd):
    # Create the PDF canvas
    pdf_canvas = canvas.Canvas(output_filename, pagesize=A4)

    # Call functions for each module you want to add
    set_cover(pdf_canvas)
    create_pdf_content(explanation_upd, visualization_upd, pdf_canvas)
    # Finalize and save the PDF
    pdf_canvas.save()

if __name__ == '__main__':
    output_filename = sys.argv[1]
    assembly = sys.argv[2]
    explanation_ori = sys.argv[3]
    visualization_ori = sys.argv[4]
    resultspath = sys.argv[5]
    group1VSgroup2 = sys.argv[6]
    
    # generate content_user.csv
    replace_dict_1 = {"sampleprefix": assembly}
    explanation_upd = explanation_ori + ".updated.csv"
    update_csv_content(explanation_ori, explanation_upd, replace_dict_1)

    # generate visualization_example_path_user.csv
    replace_dict_2 = {
        "resultspath": resultspath,
        "sampleprefix": assembly,
        "group1VSgroup2": group1VSgroup2
    }
    visualization_upd = visualization_ori + ".updated.csv"
    update_csv_content(visualization_ori,visualization_upd, replace_dict_2)

    # Constants
    page_width, page_height = A4
    HEADER_X = 20
    HEADER_Y = page_height - 20
    FOOTER_RIGHT_MARGIN = 50
    FOOTER_X = page_width - FOOTER_RIGHT_MARGIN  # 修改位置，使其总是靠右
    FOOTER_Y = 30
    FOOTER_RIGHT_MARGIN = 50
    HEADER_FONT_SIZE = 12
    FOOTER_FONT_SIZE = 10
    TEXT_FONT_SIZE = 12

    current_page_number = 1
    
    
    # Call the function to generate PDF
    generate_pdf(output_filename, explanation_upd, visualization_upd)

