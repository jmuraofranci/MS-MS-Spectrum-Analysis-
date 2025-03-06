import sys
import gzip
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from array import array
from base64 import b64decode

# Amino acid residual weights
AA_RESIDUAL_WEIGHTS = {
    'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
    'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
    'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
    'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06
}


def parse_mzxml_file(mzxml_gzip_filename, scan_number):
    """Parses an mzXML file to extract m/z values and intensity values for a given scan number."""
    expected_mz_values, intensity_values = [], []
    ns = ''
    
    with gzip.open(mzxml_gzip_filename, 'rt') as mzxml_file:
        for event, ele in ET.iterparse(mzxml_file, events=("end",)):
            if not ns:
                p = ele.tag.find('}')
                if p >= 0:
                    ns = ele.tag[:p + 1]
            
            if event == 'end' and ele.tag == ns + 'scan':
                if ele.attrib.get("num") == scan_number:
                    peaks_element = ele.find(ns + 'peaks')
                    if peaks_element is not None and peaks_element.text:
                        peaks = array('f', b64decode(peaks_element.text))
                        if sys.byteorder != 'big':
                            peaks.byteswap()
                        expected_mz_values = peaks[::2]
                        intensity_values = peaks[1::2]
                    break
                ele.clear()
    
    if not expected_mz_values:
        print("Error: Scan number does not exist in the mzXML file.")
        sys.exit(1)
    
    return expected_mz_values, intensity_values


def compute_masses(peptide_sequence):
    """Computes b-ion and y-ion masses for a given peptide sequence."""
    try:
        b_mass = [1 + sum(AA_RESIDUAL_WEIGHTS[peptide_sequence[position]] for position in range(b_ion))
                  for b_ion in range(1, len(peptide_sequence) + 1)]
        y_mass = [19 + sum(AA_RESIDUAL_WEIGHTS[peptide_sequence[position]] 
                           for position in range(len(peptide_sequence) - y_ion, len(peptide_sequence)))
                  for y_ion in range(1, len(peptide_sequence) + 1)]
        return b_mass, y_mass
    except KeyError:
        print("Error: Invalid amino acid in peptide sequence.")
        sys.exit(1)


def closest_peak(expected_mz_values, intensity_values, theoretical_mz, tolerance):
    """Finds the closest peak to a theoretical m/z value within a given tolerance."""
    intensity_match_threshold = 0.05 * max(intensity_values)
    best_peak_position = None
    best_peak_intensity = 0

    for position, mz in enumerate(expected_mz_values):
        mz_difference = abs(mz - theoretical_mz)
        if mz_difference <= tolerance and intensity_values[position] >= intensity_match_threshold:
            if best_peak_position is None or intensity_values[position] > best_peak_intensity:
                best_peak_position = position
                best_peak_intensity = intensity_values[position]
    
    return best_peak_position


def plot_spectrum(expected_mz_values, intensity_values, b_mass, y_mass, tolerance):
    """Plots the MS/MS spectrum with labeled b-ion and y-ion peaks."""
    plt.figure(figsize=(12, 8))
    plt.stem(expected_mz_values, intensity_values, linefmt='gray', markerfmt=' ', basefmt=' ')
    
    ion_dict = {'b-ions': ('blue', 'b'), 'y-ions': ('red', 'y')}
    
    for ion_list, ion_type in [(b_mass, 'b-ions'), (y_mass, 'y-ions')]:
        color, ion_label = ion_dict[ion_type]
        for position, theoretical_mz in enumerate(ion_list):
            peak_position = closest_peak(expected_mz_values, intensity_values, theoretical_mz, tolerance)
            if peak_position is not None:
                peak_mz = expected_mz_values[peak_position]
                peak_intensity = intensity_values[peak_position]
                label_position = peak_intensity + (max(intensity_values) * 0.02) if ion_type == 'b-ions' else peak_intensity
                
                plt.stem([peak_mz], [peak_intensity], linefmt=color, markerfmt=' ', basefmt=' ')
                plt.annotate(f"{ion_label}{position + 1}", xy=(peak_mz, label_position), 
                             color=color, fontsize=10, ha='center')
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    x_max = max(expected_mz_values) * 1.05
    y_max = max(intensity_values) * 1.05
    plt.xlim(0, x_max)
    plt.ylim(0, y_max)
    
    plt.xlabel('m/z')
    plt.ylabel('Relative Abundance (%)')
    plt.title('MS/MS Spectrum')
    
    plt.plot([], [], color='gray', label='Unmatched Peaks')
    plt.plot([], [], color='blue', label='b-ions Matched')
    plt.plot([], [], color='red', label='y-ions Matched')
    plt.legend(loc='upper right', fontsize=12)
    plt.show()


# Command-Line Execution
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py <mzxml_gzip_file> <scan_number> <peptide_sequence>")
        sys.exit(1)
    
    mzxml_gzip_filename = sys.argv[1]
    scan_number = sys.argv[2]
    peptide_sequence = sys.argv[3].upper()
    
    expected_mz_values, intensity_values = parse_mzxml_file(mzxml_gzip_filename, scan_number)
    b_mass, y_mass = compute_masses(peptide_sequence)
    plot_spectrum(expected_mz_values, intensity_values, b_mass, y_mass, tolerance=0.03)
