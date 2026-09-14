import os
import shutil
import SimpleITK as sitk
import nibabel as nib
import numpy as np
from totalsegmentator.python_api import totalsegmentator
from scipy.ndimage import binary_erosion, label
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd

dicom_root = r"..."
noise_folder = r"..."
results_folder = r"..."

start_index = 0

def ball_struct(radius):
    L = np.arange(-radius, radius + 1)
    X, Y, Z = np.meshgrid(L, L, L)
    return (X ** 2 + Y ** 2 + Z ** 2) <= radius ** 2

def calculate_noise_metrics(values):
    values = np.array(values)
    n_vox = len(values)
    mean_hu = np.mean(values) if n_vox > 0 else np.nan
    sd_hu = np.std(values, ddof=1) if n_vox > 1 else np.nan
    snr = mean_hu / sd_hu if sd_hu != 0 else np.nan
    return mean_hu, sd_hu, snr

def process_patient(patient_idx, patient_folder, noise_folder, results_folder, device="cpu"):

    output_path = os.path.join(patient_folder, "MASK")
    os.makedirs(output_path, exist_ok=True)

    chosen_ser = None
    chosen_series_id = None
    max_dicoms = 0
    

    for root, _, _ in os.walk(patient_folder):
        try:
            reader_tmp = sitk.ImageSeriesReader()
            series_ids_tmp = reader_tmp.GetGDCMSeriesIDs(root)

            if not series_ids_tmp:
                continue

            for series_id_tmp in series_ids_tmp:
                series_files_tmp = reader_tmp.GetGDCMSeriesFileNames(
                    root, series_id_tmp
                )
                num_dicoms = len(series_files_tmp)

                if num_dicoms > max_dicoms:
                    max_dicoms = num_dicoms
                    chosen_ser = root
                    chosen_series_id = series_id_tmp
        except Exception as e:
            print(f"Error reading {root}: {e}")
            continue

    if not chosen_ser:
        print(f"Patient {patient_idx}: DICOM-series not found")
        return

    print(f"Patient {patient_idx}: chosen {chosen_ser} ({max_dicoms} DICOM-images)")

    # DICOM to NIfTI
    nifti_output = os.path.join(output_path, f"patient{patient_idx}.nii.gz")
    reader = sitk.ImageSeriesReader()
    series_ids = reader.GetGDCMSeriesIDs(chosen_ser)
    if not series_ids:
        print(f"No serie in folder {chosen_ser}, skip")
        return
    series_files = reader.GetGDCMSeriesFileNames(chosen_ser, chosen_series_id)
    reader.SetFileNames(series_files)
    image = reader.Execute()
    sitk.WriteImage(image, nifti_output)
    print(f"NIfTI saved: {nifti_output}")

    # Read metadata directly from the DICOM images
    first_dicom = sitk.ReadImage(series_files[0])

    def dicom_value(tag, default="Unknown"):
        try:
            if first_dicom.HasMetaDataKey(tag):
                value = first_dicom.GetMetaData(tag)
                if value is not None and str(value).strip():
                    return str(value).strip()
        except Exception:
            pass
        return default

    def dicom_float(tag):
        value = dicom_value(tag, "")
        try:
            return float(value)
        except Exception:
            return np.nan

    slice_thickness = dicom_float("0018|0050")

    spacing = image.GetSpacing()
    voxel_volume = spacing[0] * spacing[1] * spacing[2]
    print(f"Voxel spacing (mm): {spacing}, Voxel volume (mm^3): {voxel_volume}")

    # TotalSegmentator
    patient_output = os.path.join(output_path, f"patient{patient_idx}_masks")
    os.makedirs(patient_output, exist_ok=True)
    try:
        totalsegmentator(
            input=nifti_output,
            output=patient_output,
            task="heartchambers_highres",
            output_type="nifti",
            quiet=False,
            device=device
        )

        print(f"Heart segmentation ready: {patient_output}")

    except Exception as e:
        print(
            f"ERROR - Patient {patient_idx}: "
            f"{os.path.basename(patient_folder)} - "
            f"TotalSegmentator segmentation failed: {e}"
        )
        return False

    # Remove all other output except the selected masks
    structures_to_keep = [
        "aorta.nii.gz",
        "heart_myocardium.nii.gz",
        "heart_ventricle_left.nii.gz",
    ]

    for f in os.listdir(patient_output):
        file_path = os.path.join(patient_output, f)
        if f not in structures_to_keep:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

    print(f"Only selected masks kept in {patient_output}")


    orig_data = sitk.GetArrayFromImage(image).transpose(1,2,0)
    
    # Myocardium
    myocardium_mask = nib.load(os.path.join(patient_output, "heart_myocardium.nii.gz")).get_fdata()
    myocardium_mask = np.rot90(myocardium_mask, k=1, axes=(0,1))
    myocardium_mask = np.flip(myocardium_mask, axis=0) > 0

    # Left ventricle mask
    lv_mask = nib.load(os.path.join(patient_output, "heart_ventricle_left.nii.gz")).get_fdata()
    lv_mask = np.rot90(lv_mask, k=1, axes=(0,1))
    lv_mask = np.flip(lv_mask, axis=0) > 0

    # Threshold HU [-200, -40] inside myocardium
    result_mask_3d = (orig_data >= -200) & (orig_data <= -40) & myocardium_mask

    # Create the outer edge mask (myocardium + LV)
    combined_mask = myocardium_mask | lv_mask
    structure = ball_struct(radius=3)

    eroded = binary_erosion(combined_mask, structure=structure)
    outer_edge = combined_mask & (~eroded)

    # Filter particles
    labeled, num_features = label(result_mask_3d)

    # Find particle labels that touch the outer edge
    edge_labels = np.unique(labeled[outer_edge])

    # Create a mask containing all other particles
    all_labels = np.arange(1, num_features+1)
    keep_labels = np.setdiff1d(all_labels, edge_labels)

    mask_clean = np.isin(labeled, keep_labels)
    print(f"Fat particles filtered")
    # Final mask: myocardium excluding the left ventricle
    result_mask_3d = mask_clean & myocardium_mask & (~lv_mask)

    # Aorta images and visualization data
    aorta_path = os.path.join(patient_output, "aorta.nii.gz")

    aorta_data = nib.load(aorta_path).get_fdata()
    aorta_data_rot = np.rot90(aorta_data, k=1, axes=(0,1))
    aorta_data_rot = np.flip(aorta_data_rot, axis=0)
    structure_aorta = ball_struct(radius=3) 
    aorta_mask_eroded = binary_erosion(aorta_data_rot, structure=structure_aorta)
    aorta_mask = aorta_mask_eroded > 0

    print(f"Start save images")

    # ==== Load all anatomical masks for visualization ====
    mask_names = [
        ("Aorta", "aorta.nii.gz"),
        ("Myocardium", "heart_myocardium.nii.gz"),
        ("Left ventricle", "heart_ventricle_left.nii.gz")
    ]

    all_masks = {}
    for mask_name, mask_file in mask_names:
        mask_path = os.path.join(patient_output, mask_file)
        if os.path.exists(mask_path):
            mask_data = nib.load(mask_path).get_fdata()
            mask_data = np.rot90(mask_data, k=1, axes=(0, 1))
            mask_data = np.flip(mask_data, axis=0) > 0
            all_masks[mask_name] = mask_data


    # Save two images for each slice
    # 1) Original CT + masks from myocardium, aorta, left ventricle, myocardium outer edge
    num_slices = orig_data.shape[2]

    mask_colors = {
        "Aorta": "blue",
        "Myocardium": "red",
        "Left ventricle": "orange",
    }

    edge_color = "green"
    fat_color = "yellow"

    for i in range(num_slices):

        orig_slice = orig_data[:, :, i]

        fig, axes = plt.subplots(1, 2, figsize=(12, 6))

        axes[0].imshow(orig_slice, cmap="gray")


        for mask_name, _ in mask_names:

            if mask_name not in mask_colors:
                continue

            if mask_name in all_masks:

                mask_slice = all_masks[mask_name][:, :, i]

                axes[0].imshow(
                    np.ma.masked_where(~mask_slice, mask_slice),
                    cmap=ListedColormap([mask_colors[mask_name]]),
                    alpha=0.35,
                    interpolation="none"
                )

        edge_slice = outer_edge[:, :, i]

        axes[0].imshow(
            np.ma.masked_where(~edge_slice, edge_slice),
            cmap=ListedColormap([edge_color]),
            alpha=0.70,
            interpolation="none"
        )

        axes[0].set_title(
            "Original CT + Myocardium + Aorta + LV + Edge"
        )
        axes[0].axis("off")

    # 2) Original CT + filtered fat particles

        axes[1].imshow(orig_slice, cmap="gray")

        particle_slice = result_mask_3d[:, :, i]

        axes[1].imshow(
            np.ma.masked_where(~particle_slice, particle_slice),
            cmap=ListedColormap([fat_color]),
            alpha=0.70,
            interpolation="none"
        )

        axes[1].set_title("Original CT + filtered fat particles")
        axes[1].axis("off")


        plt.tight_layout()

        filename = os.path.join(
            patient_output,
            f"overlay_slice_{i:03d}.png"
        )

        plt.savefig(
            filename,
            dpi=200,
            bbox_inches="tight",
            pad_inches=0.05
        )

        plt.close(fig)

    print(
        f"All {num_slices} two-panel image slices saved to folder: "
        f"{patient_output}"
    )


    # Noise assesment from the aorta region
    noise_df = pd.DataFrame(columns=[
        "Patient",
        "Region",
        "SNR",
        "Min_HU",
        "Max_HU",
        "Mean_HU",
        "SD"
    ])

    if np.any(aorta_mask):
        values = orig_data[aorta_mask]
        mean_hu, sd_hu, snr  = calculate_noise_metrics(values)

        noise_df = pd.DataFrame([{
            "Patient": patient_idx,
            "Region": "Aorta",
            "SNR": snr,
            "Min_HU": np.min(values),
            "Max_HU": np.max(values),
            "Mean_HU": mean_hu,
            "SD": sd_hu
        }])

    os.makedirs(noise_folder, exist_ok=True)
    noise_excel_path = os.path.join(
        noise_folder,
        f"Noise_control_patient_{patient_idx}.xlsx"
    )
    noise_df.to_excel(noise_excel_path, index=False)
    print(f"Aorta noise metrics saved: {noise_excel_path}")

    # Particle analysis + metal flag for possible metal artefacts
    labeled, num_features = label(result_mask_3d)

    metal_mask = orig_data >= 3000

    particle_rows = []

    for part_idx in range(1, num_features + 1):
        single_particle = (labeled == part_idx)
        num_vox = int(np.sum(single_particle))

        particle_volume = num_vox * voxel_volume
        particle_mean_hu = float(orig_data[single_particle].mean())

        # True if at least one voxel >=3000 HU occurs on any particle slice.
        slices_with_particle = np.any(single_particle, axis=(0, 1))
        contains_metal = False

        for slice_idx in np.where(slices_with_particle)[0]:
            if np.any(metal_mask[:, :, slice_idx]):
                contains_metal = True
                break

        particle_rows.append({
            "Id": len(particle_rows) + 1,
            "Number of voxels": num_vox,
            "Volume (mm3)": particle_volume,
            "Mean HU": particle_mean_hu,
            "Metal flag": bool(contains_metal)
        })

    particle_df = pd.DataFrame(
        particle_rows,
        columns=[
            "Id",
            "Number of voxels",
            "Volume (mm3)",
            "Mean HU",
            "Metal flag"
        ]
    )

    # Write results to Excel
    os.makedirs(results_folder, exist_ok=True)

    excel_path = os.path.join(
        results_folder,
        f"Results_patient_{patient_idx}.xlsx"
    )

    metadata_df = pd.DataFrame({
        "Parameter": [
            "Slice thickness",
            "Spacing",
            "Voxel volume"
        ],
        "Value": [
            slice_thickness,
            f"{spacing[0]} x {spacing[1]} x {spacing[2]} mm",
            voxel_volume
        ]
    })

    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        metadata_df.to_excel(
            writer,
            sheet_name="Results",
            index=False,
            startrow=0
        )

        particle_start_row = len(metadata_df) + 3

        particle_df.to_excel(
            writer,
            sheet_name="Results",
            index=False,
            startrow=particle_start_row
        )

    print(f"Results Excel saved: {excel_path}")



if __name__ == "__main__":

    # Sort patient folders.
    pat_folders = sorted(
        [
            f for f in os.listdir(dicom_root)
            if os.path.isdir(os.path.join(dicom_root, f))
        ]
    )

    print("=== DICOM folders ===")
    for i, folder in enumerate(pat_folders):
        print(f"Patient {i}: {folder}")

    for patient_idx, pat_folder in enumerate(
        pat_folders[start_index:], start=start_index
    ):

        patient_folder = os.path.join(dicom_root, pat_folder)

        print(f"\nProcessing {patient_idx}: {pat_folder}")

        output_path = os.path.join(patient_folder, "MASK")
        os.makedirs(output_path, exist_ok=True)

        process_patient(
            patient_idx,
            patient_folder,
            noise_folder,
            results_folder,
            device="cpu"
        )

    print("\n=== Analysis complete ===")
