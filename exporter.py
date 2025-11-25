import os 
import numpy as np

def write_results_to_txt(history_sw, history_p,grid,times,output_dir='results'):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_sw=os.path.join(output_dir,"SATURATION.txt")
    file_p=os.path.join(output_dir,"PRESSURE.txt")
    open(file_sw,'w').close()
    open(file_p,'w').close()
    print(f"Writing results to {file_sw} and {file_p}...")
    for step_idx ,(sw_matrix,p_matrix) in enumerate(zip(history_sw,history_p)):
        header = f"\n{'='*50}\nTIME STEP {step_idx} (Time: {times[step_idx]} days)\n{'='*50}\n"
        _append_matrix_to_file(file_sw,sw_matrix,grid,header,"water saturation")
        _append_matrix_to_file(file_p,p_matrix,grid,header,"Pressure (PSI)")

def _append_matrix_to_file(filename, matrix_3d, grid, header, title):
    with open(filename, 'a') as f:
        f.write(header)
        f.write(f"Property: {title}\n")
        for z in range(grid.NZ):
            f.write(f"\n layer = {z+1} Top to bottom\n")
            f.write("      ") 
            for i in range(grid.NX):
                f.write(f"{i+1:^9}") 
            f.write("\n")
            f.write("-" * (8 + 9 * grid.NX) + "\n")

            for j in range(grid.NY):
                f.write(f"J={j+1:2d} |") 
                
                for i in range(grid.NX):
                    if grid.ACT[z][j][i] == -1:
                        val_str = "    *    "
                    else:
                        val = matrix_3d[z, j, i]
                        val_str = f"{val:8.4f}" # 
                    
                    f.write(f"{val_str} ")
                f.write("\n")