"""Python script to use FS with wilson/flavio. In order to use it,
   pass the file path of the json file on the command line."""
import sys
import os.path
from wilson import Wilson
import flavio

def main():
    """Main program"""
    # Pass the path of the json file per command line
    json_file_path = sys.argv[1]
    if not os.path.exists(json_file_path):
        raise Exception('error: json file not found')
    # read Wilson coefficients from json file
    with open(json_file_path, 'r') as json_file:
        coeffs = Wilson.load_wc(json_file)
    BRXsGamma = flavio.np_prediction('BR(B->Xsgamma)', coeffs)
    print("BR(B->Xsgamma) = ", BRXsGamma)

if __name__ == "__main__":
    main()
