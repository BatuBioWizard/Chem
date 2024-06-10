import argparse
import os
from Chem_database import *
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors, SDMolSupplier, Draw
from rdkit import Chem
import tkinter as tk
from tkinter import ttk, filedialog
import csv

# database name format checker
def db_format(filename):
    if not filename.endswith('.db'):
        raise argparse.ArgumentTypeError("Your Database name must have a .db extension at the end")
    return filename


def fused_aromatic_rings_count(mol):
    # Aromatic rings (both fused and non fused)
    aromatic_rings = [ring for ring in Chem.GetSymmSSSR(mol) if mol.GetAtomWithIdx(ring[0]).GetIsAromatic()]

    # Count for fused aromatic rings
    fused_aromatic_count = 0

    # Checking for the shared atoms between aromatic rings to determine if they are fused
    for i, ring1 in enumerate(aromatic_rings):
        for ring2 in aromatic_rings[i + 1:]:
            if len(set(ring1).intersection(ring2)) > 0:
                fused_aromatic_count += 1

    return fused_aromatic_count

# Parsing SDF file
def parsing_sdf(file_path):

    suppl = SDMolSupplier(file_path)
    # Creating List to store data
    molecules_data = []

    for idx, mol in enumerate(suppl, 1):
        if mol is not None:
            # Extracting data from molecule properties
            mol_data = {
                'Molecule_Image': '',
                'Name': mol.GetProp('Name'),
                'MolWeight': float(mol.GetProp('Mol Weight')),
                'Formula': mol.GetProp('Formula'),
                'Mol_ID': int(mol.GetProp('Mol_ID')),
                'LogD': float(mol.GetProp('LogD')),
                'LogP': float(Crippen.MolLogP(mol)),
                'Donor_Count': rdMolDescriptors.CalcNumHBD(mol),
                'Acceptor_Count': rdMolDescriptors.CalcNumHBA(mol),
                'Ring_count': rdMolDescriptors.CalcNumRings(mol),
                'Aromatic_ring_count': rdMolDescriptors.CalcNumAromaticRings(mol),
                'Fused_aromatic_count': fused_aromatic_rings_count(mol),
                'PSA': rdMolDescriptors.CalcTPSA(mol),
                'Refractivity': Descriptors.MolMR(mol),
                'Rotatable_bond_counts': Descriptors.NumRotatableBonds(mol),
                'Bioavailability': 'placeholder',
                'Smiles': Chem.MolToSmiles(mol),
            }

            conditions_met = sum([
                mol_data['MolWeight'] <= 500,
                mol_data['LogP'] <= 5,
                mol_data['Donor_Count'] <= 5,
                mol_data['Acceptor_Count'] <= 10,
                mol_data['Rotatable_bond_counts'] <= 10,
                mol_data['PSA'] <= 200,
                mol_data['Fused_aromatic_count'] <= 5
            ])

            # Assign bioavailability score
            mol_data['Bioavailability'] = conditions_met
            # Appending mol_data into list
            molecules_data.append(mol_data)
            # Creating png files to store to database
            output_directory = "Molecule_Images"
            if not os.path.exists(output_directory):
                os.mkdir(output_directory)
            output_directory = "Molecule_Images"
            if not os.path.exists(output_directory):
                os.mkdir(output_directory)
            for molecule in molecules_data:
                smiles = molecule.get('Smiles')
                if smiles:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        # Define a file path for the image
                        image_file_path = f"{output_directory}/{molecule['Mol_ID']}.png"
                        img = Draw.MolToImage(mol)
                        img.save(image_file_path)
                        # Save the image file path in the molecule data
                        molecule['Molecule_Image'] = image_file_path
                    else:
                        print(f"Error: Invalid SMILES string for molecule {molecule['Mol_ID']}.")
                else:
                    print(f"Error: 'Smiles' attribute not found for molecule {molecule['Mol_ID']}.")

        else:
            print(f"Skipping molecule {idx} due to None mol object.")

    return molecules_data

class DataViewerApp:

    def __init__(self, root, data):

        self.root = root
        self.root.title("Molecule-Compound Data Viewer")
        self.data = data

        filter_title_label = ttk.Label(self.root, text="Molecule-Compound Data Viewer", font=('Arial', 40, 'bold'),
                                       foreground='white')
        filter_title_label.pack(pady=10, side='top')

        self.export_button = ttk.Button(self.root, text="Export", command=self.export_filtered_data)
        self.export_button.pack()


        # Creating Filter panel and buttons
        self.create_filter_panel()
        self.create_treeview()
        self.create_rule_buttons()

    def export_filtered_data(self):
        # Export the filtered data currently displayed in the Treeview
        filtered_data = []
        for item_id in self.tree.get_children():
            values = self.tree.item(item_id, 'values')
            filtered_data.append(values)

        # Asks user for the file name and location
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if not file_path:
            return  # User canceled the operation or closed the dialog

        # Write the filtered data to a CSV file
        with open(file_path, 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            # Write the header
            writer.writerow([column for column in self.tree['columns']])
            # Write the data
            writer.writerows(filtered_data)

    def create_filter_panel(self):

        left_frame = tk.Frame(self.root)
        left_frame.pack(side=tk.LEFT, fill=tk.Y)  # Fill the entire vertical space

        # Create the filter panel inside the left frame
        self.filter_frame = ttk.Frame(left_frame)
        self.filter_frame.pack(side=tk.LEFT, padx=10, pady=10)

        # Create a label for the filter options title within the filter frame
        filter_title_label = ttk.Label(self.filter_frame, text="Filter Options", font=('Arial', 30, 'bold'),
                                       foreground='white')
        filter_title_label.grid(row=0, column=0)

        ttk.Label(self.filter_frame, text="Select Parameter:").grid(row=1, column=0)
        self.parameter_combobox = ttk.Combobox(self.filter_frame,
                                               values=["MolWeight", "LogP", "LogD", "PSA", "Donor_Count",
                                                       "Acceptor_Count", "Ring_count", "Aromatic_ring_count",
                                                       "Fused_aromatic_count", "Refractivity", "Rotatable_bond_counts","Bioavailability"])
        self.parameter_combobox.set("Select One")  # Default parameter
        self.parameter_combobox.grid(row=1, column=1)

        self.sort_order_combobox = ttk.Combobox(self.filter_frame, values=["Asc", "Desc"])
        self.sort_order_combobox.set("Select One")  # Default sorting order
        self.sort_order_combobox.grid(row=1, column=3)

        # Max Min labels for filter data
        tk.Label(self.filter_frame, text="Minimum", font=('Helvetica', 18, 'bold')).grid(row=4, column=1)
        tk.Label(self.filter_frame, text="Maximum", font=('Helvetica', 18, 'bold')).grid(row=4, column=3)

        # Creating Filter objects that we want to use
        ttk.Label(self.filter_frame, text="MolWeight:").grid(row=5, column=0)
        self.molweight_min_entry = ttk.Entry(self.filter_frame, width=8)
        self.molweight_min_entry.grid(row=5, column=1)
        ttk.Label(self.filter_frame, text="to").grid(row=5, column=2)
        self.molweight_max_entry = ttk.Entry(self.filter_frame, width=8)
        self.molweight_max_entry.grid(row=5, column=3)

        ttk.Label(self.filter_frame, text="logP:").grid(row=6, column=0)
        self.logp_min_entry = ttk.Entry(self.filter_frame, width=8)
        self.logp_min_entry.grid(row=6, column=1)
        ttk.Label(self.filter_frame, text="to").grid(row=6, column=2)
        self.logp_max_entry = ttk.Entry(self.filter_frame, width=8)
        self.logp_max_entry.grid(row=6, column=3)

        ttk.Label(self.filter_frame, text="logD:").grid(row=7, column=0)
        self.logd_min_entry = ttk.Entry(self.filter_frame, width=8)
        self.logd_min_entry.grid(row=7, column=1)
        ttk.Label(self.filter_frame, text="to").grid(row=7, column=2)
        self.logd_max_entry = ttk.Entry(self.filter_frame, width=8)
        self.logd_max_entry.grid(row=7, column=3)

        ttk.Label(self.filter_frame, text="PSA:").grid(row=8, column=0)
        self.psa_min_entry = ttk.Entry(self.filter_frame, width=8)
        self.psa_min_entry.grid(row=8, column=1)
        ttk.Label(self.filter_frame, text="to").grid(row=8, column=2)
        self.psa_max_entry = ttk.Entry(self.filter_frame, width=8)
        self.psa_max_entry.grid(row=8, column=3)

        ttk.Label(self.filter_frame, text="Donor_Count:").grid(row=9, column=0)
        self.donor_min_entry = ttk.Entry(self.filter_frame, width=8)
        self.donor_min_entry.grid(row=9, column=1)
        ttk.Label(self.filter_frame, text="to").grid(row=9, column=2)
        self.donor_max_entry = ttk.Entry(self.filter_frame, width=8)
        self.donor_max_entry.grid(row=9, column=3)

        ttk.Label(self.filter_frame, text="Acceptor_Count:").grid(row=10, column=0)
        self.acceptor_min_entry = ttk.Entry(self.filter_frame, width=8)
        self.acceptor_min_entry.grid(row=10, column=1)
        ttk.Label(self.filter_frame, text="to").grid(row=10, column=2)
        self.acceptor_max_entry = ttk.Entry(self.filter_frame, width=8)
        self.acceptor_max_entry.grid(row=10, column=3)

        ttk.Label(self.filter_frame, text="Ring_count:").grid(row=11, column=0)
        self.ring_min_entry = ttk.Entry(self.filter_frame, width=8)
        self.ring_min_entry.grid(row=11, column=1)
        ttk.Label(self.filter_frame, text="to").grid(row=11, column=2)
        self.ring_max_entry = ttk.Entry(self.filter_frame, width=8)
        self.ring_max_entry.grid(row=11, column=3)

        ttk.Label(self.filter_frame, text="Aromatic_ring_count:").grid(row=12, column=0)
        self.aromatic_ring_min_entry = ttk.Entry(self.filter_frame, width=8)
        self.aromatic_ring_min_entry.grid(row=12, column=1)
        ttk.Label(self.filter_frame, text="to").grid(row=12, column=2)
        self.aromatic_ring_max_entry = ttk.Entry(self.filter_frame, width=8)
        self.aromatic_ring_max_entry.grid(row=12, column=3)

        ttk.Label(self.filter_frame, text="Fused_aromatic_count:").grid(row=13, column=0)
        self.fused_aromatic_min_entry = ttk.Entry(self.filter_frame, width=8)
        self.fused_aromatic_min_entry.grid(row=13, column=1)
        ttk.Label(self.filter_frame, text="to").grid(row=13, column=2)
        self.fused_aromatic_max_entry = ttk.Entry(self.filter_frame, width=8)
        self.fused_aromatic_max_entry.grid(row=13, column=3)

        ttk.Label(self.filter_frame, text="Refractivity:").grid(row=14, column=0)
        self.refractivity_min_entry = ttk.Entry(self.filter_frame, width=8)
        self.refractivity_min_entry.grid(row=14, column=1)
        ttk.Label(self.filter_frame, text="to").grid(row=14, column=2)
        self.refractivity_max_entry = ttk.Entry(self.filter_frame, width=8)
        self.refractivity_max_entry.grid(row=14, column=3)

        ttk.Label(self.filter_frame, text="Rotatable_bond_counts:").grid(row=15, column=0)
        self.rotatable_bond_min_entry = ttk.Entry(self.filter_frame, width=8)
        self.rotatable_bond_min_entry.grid(row=15, column=1)
        ttk.Label(self.filter_frame, text="to").grid(row=15, column=2)
        self.rotatable_bond_max_entry = ttk.Entry(self.filter_frame, width=8)
        self.rotatable_bond_max_entry.grid(row=15, column=3)

        # Adding Filter button
        self.filter_button = ttk.Button(self.filter_frame, text="Filter", command=self.filter_data)
        self.filter_button.grid(row=16, column=1)
        # Adding Clear button
        self.clear_button = ttk.Button(self.filter_frame, text="Clear", command=self.clear_filter)
        self.clear_button.grid(row=16, column=3)
        # Adding Filter Count
        ttk.Label(self.filter_frame, text="Filtered Count:", font=('Helvetica', 18, 'bold')).grid(row=17, column=0)
        self.filtered_count_entry = ttk.Entry(self.filter_frame, width=8)
        self.filtered_count_entry.grid(row=17, column=1)
        # Adding Lipinski's Rule Button
        self.lipinski_button = ttk.Button(self.filter_frame, text="Apply Lipinski's Rule",
                                          command=self.lipinski_rule, width=20)
        self.lipinski_button.grid(row=18, column=0, padx=5, pady=5)
        # Adding Lead-likeness Rule Button
        self.leadlikeness_button = ttk.Button(self.filter_frame, text="Apply Lead-likeness Rule",
                                              command=self.lead_likeness_rule, width=20)
        self.leadlikeness_button.grid(row=19, column=0, padx=5, pady=5)
        # Adding Bioavailability Rule Button
        self.bioavailability_button = ttk.Button(self.filter_frame, text="Apply Bioavailability Rule",
                                                 command=self.bioavailability_rule, width=20)
        self.bioavailability_button.grid(row=20, column=0, padx=5, pady=5)

    def clear_filter(self):
        # Clearing all entry fields
        self.molweight_min_entry.delete(0, tk.END)
        self.molweight_max_entry.delete(0, tk.END)
        self.logp_min_entry.delete(0, tk.END)
        self.logp_max_entry.delete(0, tk.END)
        self.logd_min_entry.delete(0, tk.END)
        self.logd_max_entry.delete(0, tk.END)
        self.psa_min_entry.delete(0, tk.END)
        self.psa_max_entry.delete(0, tk.END)
        self.donor_min_entry.delete(0, tk.END)
        self.donor_max_entry.delete(0, tk.END)
        self.acceptor_min_entry.delete(0, tk.END)
        self.acceptor_max_entry.delete(0, tk.END)
        self.ring_min_entry.delete(0, tk.END)
        self.ring_max_entry.delete(0, tk.END)
        self.aromatic_ring_min_entry.delete(0, tk.END)
        self.aromatic_ring_max_entry.delete(0, tk.END)
        self.fused_aromatic_min_entry.delete(0, tk.END)
        self.fused_aromatic_max_entry.delete(0, tk.END)
        self.refractivity_min_entry.delete(0, tk.END)
        self.refractivity_max_entry.delete(0, tk.END)
        self.rotatable_bond_min_entry.delete(0, tk.END)
        self.rotatable_bond_max_entry.delete(0, tk.END)
        self.search_entry.delete(0, tk.END)
        self.filtered_count_entry.delete(0, tk.END)
        # Destroy existing tree and reload data to have initial look
        self.tree.delete(*self.tree.get_children())
        self.load_data()

    def create_treeview(self):

        self.tree = ttk.Treeview(self.root)
        self.tree["columns"] = tuple(self.data[0].keys())

        for column in self.tree["columns"]:
            self.tree.heading(column, text=column)
            # Adjusting column widths
            if column == "Name":
                self.tree.column(column, width=600)  # Make "Name" column wider
            else:
                self.tree.column(column, width=220)  # Set other columns width

        # Create a horizontal scrollbar
        h_scrollbar = ttk.Scrollbar(self.root, orient="horizontal", command=self.tree.xview)

        # Configure the treeview to use the horizontal scrollbar
        self.tree.configure(xscrollcommand=h_scrollbar.set)

        # Pack the treeview and scrollbar
        self.tree.pack(expand=True, fill="both")
        h_scrollbar.pack(side="bottom", fill="x")

        self.load_data()


    def load_data(self):
        for item in self.data:
            self.tree.insert("", "end", values=tuple(item.values()))

    def filter_data(self):
        # Setting combo box for ranking
        selected_parameter = self.parameter_combobox.get()
        selected_sort_order = self.sort_order_combobox.get()

        # Geting the values from the entry fields
        molweight_min = float(self.molweight_min_entry.get()) if self.molweight_min_entry.get() else None
        molweight_max = float(self.molweight_max_entry.get()) if self.molweight_max_entry.get() else None
        logp_min = float(self.logp_min_entry.get()) if self.logp_min_entry.get() else None
        logp_max = float(self.logp_max_entry.get()) if self.logp_max_entry.get() else None
        logd_min = float(self.logd_min_entry.get()) if self.logd_min_entry.get() else None
        logd_max = float(self.logd_max_entry.get()) if self.logd_max_entry.get() else None
        psa_min = float(self.psa_min_entry.get()) if self.psa_min_entry.get() else None
        psa_max = float(self.psa_max_entry.get()) if self.psa_max_entry.get() else None
        donor_min = int(self.donor_min_entry.get()) if self.donor_min_entry.get() else None
        donor_max = int(self.donor_max_entry.get()) if self.donor_max_entry.get() else None
        acceptor_min = int(self.acceptor_min_entry.get()) if self.acceptor_min_entry.get() else None
        acceptor_max = int(self.acceptor_max_entry.get()) if self.acceptor_max_entry.get() else None
        ring_min = int(self.ring_min_entry.get()) if self.ring_min_entry.get() else None
        ring_max = int(self.ring_max_entry.get()) if self.ring_max_entry.get() else None
        aromatic_ring_min = int(self.aromatic_ring_min_entry.get()) if self.aromatic_ring_min_entry.get() else None
        aromatic_ring_max = int(self.aromatic_ring_max_entry.get()) if self.aromatic_ring_max_entry.get() else None
        fused_aromatic_min = int(self.fused_aromatic_min_entry.get()) if self.fused_aromatic_min_entry.get() else None
        fused_aromatic_max = int(self.fused_aromatic_max_entry.get()) if self.fused_aromatic_max_entry.get() else None
        refractivity_min = float(self.refractivity_min_entry.get()) if self.refractivity_min_entry.get() else None
        refractivity_max = float(self.refractivity_max_entry.get()) if self.refractivity_max_entry.get() else None
        rotatable_bond_min = int(self.rotatable_bond_min_entry.get()) if self.rotatable_bond_min_entry.get() else None
        rotatable_bond_max = int(self.rotatable_bond_max_entry.get()) if self.rotatable_bond_max_entry.get() else None

        # Filtering the data based on the entered values
        filtered_data = [item for item in self.data if
                         (molweight_min is None or molweight_min <= item['MolWeight']) and
                         (molweight_max is None or item['MolWeight'] <= molweight_max) and
                         (logp_min is None or logp_min <= item['LogP']) and
                         (logp_max is None or item['LogP'] <= logp_max) and
                         (logd_min is None or logd_min <= item['LogD']) and
                         (logd_max is None or item['LogD'] <= logd_max) and
                         (psa_min is None or psa_min <= item['PSA']) and
                         (psa_max is None or item['PSA'] <= psa_max) and
                         (donor_min is None or donor_min <= item['Donor_Count']) and
                         (donor_max is None or item['Donor_Count'] <= donor_max) and
                         (acceptor_min is None or acceptor_min <= item['Acceptor_Count']) and
                         (acceptor_max is None or item['Acceptor_Count'] <= acceptor_max) and
                         (ring_min is None or ring_min <= item['Ring_count']) and
                         (ring_max is None or item['Ring_count'] <= ring_max) and
                         (aromatic_ring_min is None or aromatic_ring_min <= item['Aromatic_ring_count']) and
                         (aromatic_ring_max is None or item['Aromatic_ring_count'] <= aromatic_ring_max) and
                         (fused_aromatic_min is None or fused_aromatic_min <= item['Fused_aromatic_count']) and
                         (fused_aromatic_max is None or item['Fused_aromatic_count'] <= fused_aromatic_max) and
                         (refractivity_min is None or refractivity_min <= item['Refractivity']) and
                         (refractivity_max is None or item['Refractivity'] <= refractivity_max) and
                         (rotatable_bond_min is None or rotatable_bond_min <= item['Rotatable_bond_counts']) and
                         (rotatable_bond_max is None or item['Rotatable_bond_counts'] <= rotatable_bond_max)]

        # Sort the filtered data based on the selected parameter and sort order
        if selected_sort_order == 'Asc':
            filtered_data.sort(key=lambda x: x[selected_parameter])
        elif selected_sort_order == 'Desc':
            filtered_data.sort(key=lambda x: x[selected_parameter], reverse=True)

        # Clearing the current data in the treeview
        self.tree.delete(*self.tree.get_children())

        # Inserting the filtered data into the treeview
        for item in filtered_data:
            self.tree.insert("", "end", values=tuple(item.values()))

        # Checking Filtered Count
        filtered_count = len(self.tree.get_children())
        self.filtered_count_entry.delete(0, tk.END)
        self.filtered_count_entry.insert(0, filtered_count)

    def create_rule_buttons(self):

        self.rule_button_frame = ttk.Frame(self.root)
        self.rule_button_frame.pack(pady=10)

        # Creating search bar label
        ttk.Label(self.rule_button_frame, text="Search Molecule&Smile Name:").grid(row=0, column=0, sticky="w")
        self.search_entry = ttk.Entry(self.rule_button_frame, width=30)
        self.search_entry.grid(row=0, column=1)

        # Add a search button
        self.search_button = ttk.Button(self.rule_button_frame, text="Search", command=self.search_molecule)
        self.search_button.grid(row=0, column=2)

    def display_filter_criteria(self):

        # Display filter criteria for Lipinski Rule in entry widgets

        # For MolWeight
        self.molweight_min_entry.delete(0, tk.END)
        self.molweight_min_entry.insert(0, "0")
        self.molweight_max_entry.delete(0, tk.END)
        self.molweight_max_entry.insert(0, "500")

        # For LogP
        self.logp_min_entry.delete(0, tk.END)
        self.logp_min_entry.insert(0, "")
        self.logp_max_entry.delete(0, tk.END)
        self.logp_max_entry.insert(0, "5")

        # For Donor_Count
        self.donor_min_entry.delete(0, tk.END)
        self.donor_min_entry.insert(0, "0")
        self.donor_max_entry.delete(0, tk.END)
        self.donor_max_entry.insert(0, "5")

        # For Acceptor_Count
        self.acceptor_min_entry.delete(0, tk.END)
        self.acceptor_min_entry.insert(0, "0")
        self.acceptor_max_entry.delete(0, tk.END)
        self.acceptor_max_entry.insert(0, "10")

        # Checking Filtered Count
        filtered_count = len(self.tree.get_children())
        self.filtered_count_entry.delete(0, tk.END)
        self.filtered_count_entry.insert(0, filtered_count)

    def display_filter_criteria_2(self):

        # Display filter criteria for Lead-likeness Rule in entry widgets

        # For MolWeight
        self.molweight_min_entry.delete(0, tk.END)
        self.molweight_min_entry.insert(0, "0")
        self.molweight_max_entry.delete(0, tk.END)
        self.molweight_max_entry.insert(0, "450")

        # For LogD
        self.logd_min_entry.delete(0, tk.END)
        self.logd_min_entry.insert(0, "-4")
        self.logd_max_entry.delete(0, tk.END)
        self.logd_max_entry.insert(0, "4")

        # For Ring_count
        self.ring_min_entry.delete(0, tk.END)
        self.ring_min_entry.insert(0, "0")
        self.ring_max_entry.delete(0, tk.END)
        self.ring_max_entry.insert(0, "4")

        # For Rotatable_bond_counts
        self.rotatable_bond_min_entry.delete(0, tk.END)
        self.rotatable_bond_min_entry.insert(0, "0")
        self.rotatable_bond_max_entry.delete(0, tk.END)
        self.rotatable_bond_max_entry.insert(0, "10")

        # For Donor_Count
        self.donor_min_entry.delete(0, tk.END)
        self.donor_min_entry.insert(0, "0")
        self.donor_max_entry.delete(0, tk.END)
        self.donor_max_entry.insert(0, "5")

        # For Acceptor_Count
        self.acceptor_min_entry.delete(0, tk.END)
        self.acceptor_min_entry.insert(0, "0")
        self.acceptor_max_entry.delete(0, tk.END)
        self.acceptor_max_entry.insert(0, "8")

        # Checking Filtered Count
        filtered_count = len(self.tree.get_children())
        self.filtered_count_entry.delete(0, tk.END)
        self.filtered_count_entry.insert(0, filtered_count)

    def display_filter_criteria_3(self):

        # Display filter criteria for Bioavailability Rule in entry widgets

        # For MolWeight
        self.molweight_min_entry.delete(0, tk.END)
        self.molweight_min_entry.insert(0, "0")
        self.molweight_max_entry.delete(0, tk.END)
        self.molweight_max_entry.insert(0, "500")

        # For LogP
        self.logp_min_entry.delete(0, tk.END)
        self.logp_min_entry.insert(0, "")
        self.logp_max_entry.delete(0, tk.END)
        self.logp_max_entry.insert(0, "5")

        # For Donor_Count
        self.donor_min_entry.delete(0, tk.END)
        self.donor_min_entry.insert(0, "0")
        self.donor_max_entry.delete(0, tk.END)
        self.donor_max_entry.insert(0, "5")

        # For Acceptor_Count
        self.acceptor_min_entry.delete(0, tk.END)
        self.acceptor_min_entry.insert(0, "0")
        self.acceptor_max_entry.delete(0, tk.END)
        self.acceptor_max_entry.insert(0, "10")

        # For Rotatable_bond_counts
        self.rotatable_bond_min_entry.delete(0, tk.END)
        self.rotatable_bond_min_entry.insert(0, "0")
        self.rotatable_bond_max_entry.delete(0, tk.END)
        self.rotatable_bond_max_entry.insert(0, "10")

        # For PSA
        self.psa_min_entry.delete(0, tk.END)
        self.psa_min_entry.insert(0, "0")
        self.psa_max_entry.delete(0, tk.END)
        self.psa_max_entry.insert(0, "200")

        # For Fused_aromatic_count
        self.fused_aromatic_min_entry.delete(0, tk.END)
        self.fused_aromatic_min_entry.insert(0, "0")
        self.fused_aromatic_max_entry.delete(0, tk.END)
        self.fused_aromatic_max_entry.insert(0, "5")

        # Checking Filtered Count
        filtered_count = len(self.tree.get_children())
        self.filtered_count_entry.delete(0, tk.END)
        self.filtered_count_entry.insert(0, filtered_count)

    def lipinski_rule(self):
        # Clearing entry table before searching
        self.clear_filter()
        filtered_data = [item for item in self.data if
                         item['MolWeight'] <= 500 and
                         item['LogP'] <= 5 and
                         item['Donor_Count'] <= 5 and
                         item['Acceptor_Count'] <= 10]

        # Sorting the filtered data based on MolWeight in descending order
        filtered_data.sort(key=lambda x: x['MolWeight'], reverse=True)

        # Clearing the current data in the treeview
        self.tree.delete(*self.tree.get_children())

        # Inserting the filtered data into the treeview
        for item in filtered_data:
            self.tree.insert("", "end", values=tuple(item.values()))
        # Calling Display function
        self.display_filter_criteria()


    def lead_likeness_rule(self):
        # Clearing entry table before searching
        self.clear_filter()
        filtered_data = [item for item in self.data if
                         item['MolWeight'] <= 450 and
                         -4 <= item['LogD'] <= 4 and
                         item['Ring_count'] <= 4 and
                         item['Rotatable_bond_counts'] <= 10 and
                         item['Donor_Count'] <= 5 and
                         item['Acceptor_Count'] <= 8]

        # Sorting the filtered data based on MolWeight in descending order
        filtered_data.sort(key=lambda x: x['MolWeight'], reverse=True)

        # Clearing the current data in the treeview
        self.tree.delete(*self.tree.get_children())

        # Inserting the filtered data into the treeview
        for item in filtered_data:
            self.tree.insert("", "end", values=tuple(item.values()))
        # Calling Display function
        self.display_filter_criteria_2()


    def bioavailability_rule(self):
        # Clearing entry table before searching
        self.clear_filter()
        count_conditions_met = lambda item: sum([
            item['MolWeight'] <= 500,
            item['LogP'] <= 5,
            item['Donor_Count'] <= 5,
            item['Acceptor_Count'] <= 10,
            item['Rotatable_bond_counts'] <= 10,
            item['PSA'] <= 200,
            item['Fused_aromatic_count'] <= 5
        ])

        filtered_data = [item for item in self.data if count_conditions_met(item) >= 6]

        # Sorting the filtered data based on MolWeight in descending order
        filtered_data.sort(key=lambda x: x['MolWeight'], reverse=True)

        # Clearing the current data in the treeview
        self.tree.delete(*self.tree.get_children())

        # Inserting the filtered data into the treeview
        for item in filtered_data:
            self.tree.insert("", "end", values=tuple(item.values()))
        # Calling Display function
        self.display_filter_criteria_3()

    def search_molecule(self):
        # Get the search query from the entry widget
        query = self.search_entry.get().strip().lower()

        # Clear the current data in the treeview
        for i in self.tree.get_children():
            self.tree.delete(i)

        # Filter the data to include molecules with names or 'Smiles' containing the search query
        filtered_data = [item for item in self.data if
                         query in item['Name'].lower() or
                         query in item['Smiles'].lower()]

        filtered_data.sort(key=lambda x: x['MolWeight'], reverse=True)

        # Insert the filtered data into the treeview
        for item in filtered_data:
            self.tree.insert("", "end", values=tuple(item.values()))

        # Checking Filtered Count
        filtered_count = len(self.tree.get_children())
        self.filtered_count_entry.delete(0, tk.END)
        self.filtered_count_entry.insert(0, filtered_count)
