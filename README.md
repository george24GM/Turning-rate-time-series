# **README**  

## **Repository Overview**  
This repository contains MATLAB scripts for the simulations of the paper. The scripts generate plots and numerical results used in Figures 3 and 4 and retrieve values for Table 1.  

**Note:** To ensure proper execution, all files must be placed in the same directory.  

---

## **File Descriptions**  

### **1. `Setting_B.m`**  
- Generates the plots for **Figure 4**, which depict the power of the test for detecting changes in the AR(1) parameter.  
- By printing the variable **`v`** in the MATLAB terminal, users can retrieve the values corresponding to **Table 1 (last row)**.  
- The innovations (error terms) can be modified automatically by commenting/uncommenting the relevant distribution in the script.  

### **2. `Setting_A.m`**  
- Generates the plots for **Figure 3** related to the study.  

### **3. `EEG_reader.m`**  
- Processes EEG data stored in **.edf format**, downloaded from **PhysioNet**, to detect potential change points.  

### **4. `Detect_Change_Point.m`**  
- Implements an algorithm to detect change points in a given time series.  
- The generated plots are specifically tailored to the **real data example** presented in the paper.  

---
