import subprocess
import json
import os

new_json = {}

with open("smiles_by_ref.json", "r") as f:
    json_data = json.load(f)

for count, item in enumerate(list(json_data)):
    print(item)
    if count % 100 == 0:
        print("Iterate through item {}".format(count))
    result = subprocess.run("obabel -:\"{}\" -oxyz --gen3d -O output.xyz".format(item), shell=True)
    if result.returncode != 0:
        print(f"Open Babel failed for item: {item}")
        continue

    xyz_file = "output.xyz"
    if os.path.isfile(xyz_file):
        with open(xyz_file, "r") as xyz_f:
            xyz_content = xyz_f.read()
    else:
        print(item)
    # Add XYZ content to original dictionary
    new_json[item] = {}
    new_json[item]["xyz_coordinate"] = xyz_content
    new_json[item]["ref_code"] = json_data[item]
    subprocess.run("rm {}".format(xyz_file), shell=True)
# Save updated JSON data
with open("test.json", "w") as f:
    json.dump(new_json, f, indent=4)

