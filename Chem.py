from Chem_functions import *
def main():
    # Commandline prep
    parser = argparse.ArgumentParser(description="Chem5042 Assessment 2", allow_abbrev=False)
    parser.add_argument("db_name", type=str, help="Database name")
    parser.add_argument("--createdb", action="store_true", help="Creating the database structure")
    parser.add_argument("sdfFile", help='Please enter your sdf file after --sdfFile input')
    args = parser.parse_args()
    # Creating Path for SDF file to parse
    sdf_file = args.sdfFile
    molecules_data = parsing_sdf(sdf_file)
    # Whenever --createdb used in commandline creating db and inserting data in it
    if args.createdb:
        db_creation(args.db_name)
        for molecule_data in molecules_data:
            image_path = molecule_data['Molecule_Image']
            inserting_db(args.db_name, molecule_data, image_path)
    # Creating UI
    root = tk.Tk()
    root.geometry("1000x600")
    app = DataViewerApp(root, molecules_data)
    root.mainloop()

if __name__ == "__main__":
    main()

# python3 Chem.py Chem.db --createdb Molecules1.sdf