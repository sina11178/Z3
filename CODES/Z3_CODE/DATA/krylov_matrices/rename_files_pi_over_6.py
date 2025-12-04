import os

# Path to the directory containing your files

folder = r"C:\Users\frsin\OneDrive\Desktop\Grad school stuff\Research_Notes\Kolodrubetz Group\Z3\Z3 Project\Z3\CODES\Z3_CODE\DATA\krylov_matrices\krylov_matrices_120425_vary_gamma"

# Loop through all files recursively
for root, dirs, files in os.walk(folder):
    for filename in files:
        if "0.5235987755982988" in filename:
            new_filename = filename.replace("0.5235987755982988", "0.524")
            old_path = os.path.join(root, filename)
            new_path = os.path.join(root, new_filename)
            os.rename(old_path, new_path)
            print(f"Renamed:\n{old_path}\n→ {new_path}\n")