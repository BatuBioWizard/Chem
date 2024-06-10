import sqlite3

# Database creation
def db_creation(db_name):
    database = """
    CREATE TABLE IF NOT EXISTS Molecules (
    id INTEGER PRIMARY KEY,
    Molecule_Image BLOB,  
    Name TEXT,
    MolWeight REAL,
    Formula TEXT,
    Mol_ID INTEGER,
    LogD REAL,
    LogP REAL,
    Donor_Count INTEGER,
    Acceptor_Count INTEGER,
    Ring_count INTEGER,
    Aromatic_ring_count INTEGER,
    Fused_aromatic_count INTEGER,
    PSA REAL,
    Refractivity REAL,
    Rotatable_bond_counts INTEGER,
    Bioavailability INTEGER
);
"""
    with sqlite3.connect(db_name) as db:
       cursor = db.cursor()
       cursor.executescript(database)
       db.commit()

# Converting strings to binary data
def convert_image_to_binary(image_path):
    with open(image_path, 'rb') as file:
        blob_data = file.read()
    return blob_data
# Inserting data
def inserting_db(db_name, molecule_data, image_path):
    try:
        with sqlite3.connect(db_name) as db:
            image_blob = convert_image_to_binary(image_path)
            cursor = db.cursor()
            # Checking if the molecule with same name or identifier already exists
            cursor.execute("SELECT * FROM Molecules WHERE Name = ? OR Mol_ID = ?",
                           (molecule_data['Name'], molecule_data['Mol_ID']))
            existing_record = cursor.fetchone()
            if existing_record:
                # Updating existing record
                cursor.execute(
                    """UPDATE Molecules SET Molecule_Image = ?, MolWeight = ?, Formula = ?, LogD = ?, LogP = ?,
                    Donor_Count = ?, Acceptor_Count = ?, Ring_count = ?, Aromatic_ring_count = ?,
                    Fused_aromatic_count = ?, PSA = ?, Refractivity = ?, Rotatable_bond_counts = ?, Bioavailability = ?
                    WHERE id = ?""",
                    (image_blob, molecule_data['MolWeight'], molecule_data['Formula'], molecule_data['LogD'], molecule_data['LogP'],
                    molecule_data['Donor_Count'], molecule_data['Acceptor_Count'], molecule_data['Ring_count'],
                    molecule_data['Aromatic_ring_count'], molecule_data['Fused_aromatic_count'], molecule_data['PSA'],
                    molecule_data['Refractivity'], molecule_data['Rotatable_bond_counts'], molecule_data['Bioavailability'], existing_record[0])
                )
            else:
                # Inserting new record
                cursor.execute(
                    """INSERT INTO Molecules (Molecule_Image, Name, MolWeight, Formula, Mol_ID, LogD, LogP,
                    Donor_Count, Acceptor_Count, Ring_count, Aromatic_ring_count, Fused_aromatic_count,
                    PSA, Refractivity, Rotatable_bond_counts, Bioavailability)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (image_blob, molecule_data['Name'], molecule_data['MolWeight'], molecule_data['Formula'], molecule_data['Mol_ID'],
                    molecule_data['LogD'], molecule_data['LogP'], molecule_data['Donor_Count'], molecule_data['Acceptor_Count'],
                    molecule_data['Ring_count'], molecule_data['Aromatic_ring_count'], molecule_data['Fused_aromatic_count'],
                    molecule_data['PSA'], molecule_data['Refractivity'], molecule_data['Rotatable_bond_counts'], molecule_data['Bioavailability'])
                )
            db.commit()
    except sqlite3.Error as e:
        print("Error, you failed to insert into the database:", e)


