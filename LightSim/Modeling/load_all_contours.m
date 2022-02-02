% load contours
global model_C_data_p5mm
global model_C_data_1mm
global model_C_data_1p5mm
global model_C_data_2mm

model_C_data_p5mm = load_and_clean_contours("D:\data\pt\gridsearchbest\gridsearchbestvalid_p5mm\gridsearchbestvalid_p5mm.mc2", 2);
model_C_data_1mm = load_and_clean_contours("D:\data\pt\gridsearchbest\gridsearchbestvalid_1mm\gridsearchbestvalid_1mm.mc2", 2);
model_C_data_1p5mm = load_and_clean_contours("D:\data\pt\gridsearchbest\gridsearchbestvalid_1p5mm\gridsearchbestvalid_1p5mm.mc2", 2);
model_C_data_2mm = load_and_clean_contours("D:\data\pt\gridsearchbest\gridsearchbestvalid_2mm\gridsearchbestvalid_2mm.mc2", 2);