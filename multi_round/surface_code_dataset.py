UP = lambda x: (x[0]-1, x[1])
DW = lambda x: (x[0]+1, x[1])
LF = lambda x: (x[0], x[1]-1)
RF = lambda x: (x[0], x[1]+1)
import numpy as np
import pymatching
from .surface_code import *

def disassemble_sample(test_sample, error_type):
    # z error: 0, x_error: 1
    syndrome_pos = 2*error_type
    detection_event = [k for k in test_sample[0][syndrome_pos]]
    for test_round in range(1, len(test_sample)):
        # detection event are triggered if any signal is flipped.
        detection_event += list_xor(test_sample[test_round-1][syndrome_pos], test_sample[test_round][syndrome_pos])
    noise = [k for k in test_sample[len(test_sample)-1][syndrome_pos+1]]
    return detection_event, noise

# generating dataset
def gene_mwpm_dataset(code, num_round, err_prob, num_samples, retrieve_samples=None):
    data_qubits = code.data_qubits
    mz_qubits = code.mz_qubits
    mx_qubits = code.mx_qubits
    distance = code.distance
    data = code._sample(num_round, err_prob, num_samples)
    if retrieve_samples != None:
        retrieve_samples.extend(data)
    xerr_matching = code.get_z_signal_matching_graph(rounds=num_round)
    zerr_matching = code.get_x_signal_matching_graph(rounds=num_round)
    x_error_type = 1
    z_error_type = 0
    samples = []
    labels = []
    for i in range(num_samples):
        signal_by_xerr, x_noise = disassemble_sample(data[i], x_error_type)
        signal_by_zerr, z_noise = disassemble_sample(data[i], z_error_type)
        xerr_prediction = xerr_matching.decode(signal_by_xerr) # decoding x noise
        zerr_prediction = zerr_matching.decode(signal_by_zerr) # decoding z noise
        # generating syndrome_by_x_err
        z_signal_images = []
        cnt = 0
        for k in range(num_round):
            image = []
            for j in range(0,2*distance-1):
                row = []
                for i in range(0, 2*distance-1):
                    if (j,i) in mz_qubits:
                        row.append(signal_by_xerr[cnt])
                        cnt += 1
                    else:
                        row.append(0)
                image.append(row)
            z_signal_images.append(image)
        cnt = 0
        x_signal_images = []
        for k in range(num_round):
            image = []
            for j in range(0,2*distance-1):
                row = []
                for i in range(0, 2*distance-1):
                    if (j,i) in mx_qubits:
                        row.append(signal_by_zerr[cnt])
                        cnt += 1
                    else:
                        row.append(0)
                image.append(row)
            x_signal_images.append(image)
        z_boundary_image = []
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                if j == 0 or j == 2*distance-2:
                    if i % 2 == 0:
                        row.append(1)
                    else:
                        row.append(0)
                else:
                    row.append(0)
            z_boundary_image.append(row)
        x_boundary_image = []
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                if i == 0 or i == 2*distance-2:
                    if j % 2 == 0:
                        row.append(1)
                    else:
                        row.append(0)
                else:
                    row.append(0)
            x_boundary_image.append(row)
        x_err_pred_image = []
        cnt = 0
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                if (j,i) in data_qubits:
                    row.append(xerr_prediction[cnt])
                    cnt += 1
                else:
                    row.append(0)
            x_err_pred_image.append(row)
        z_err_pred_image = []
        cnt = 0
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                if (j,i) in data_qubits:
                    row.append(zerr_prediction[cnt])
                    cnt += 1
                else:
                    row.append(0)
            z_err_pred_image.append(row)
        z_mea_image = []
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                row.append(0)
            z_mea_image.append(row)
        x_mea_image = []
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                row.append(0)
            x_mea_image.append(row)
        samples.append(z_signal_images+x_signal_images+[z_boundary_image, x_boundary_image])
        labels.append([x_err_pred_image, z_err_pred_image, z_mea_image, x_mea_image])
    return np.array(samples), np.array(labels)

def gene_dataset(code, num_round, err_prob, num_samples, retrieve_samples=None):
    data_qubits = code.data_qubits
    mz_qubits = code.mz_qubits
    mx_qubits = code.mx_qubits
    distance = code.distance
    data = code._sample(num_round, err_prob, num_samples)
    if retrieve_samples != None:
        retrieve_samples.extend(data)
    x_error_type = 1
    z_error_type = 0
    samples = []
    labels = []
    for i in range(num_samples):
        signal_by_xerr, x_noise = disassemble_sample(data[i], x_error_type)
        signal_by_zerr, z_noise = disassemble_sample(data[i], z_error_type)
        # generating syndrome_by_x_err
        z_signal_images = []
        cnt = 0
        for k in range(num_round):
            image = []
            for j in range(0,2*distance-1):
                row = []
                for i in range(0, 2*distance-1):
                    if (j,i) in mz_qubits:
                        row.append(signal_by_xerr[cnt])
                        cnt += 1
                    else:
                        row.append(0)
                image.append(row)
            z_signal_images.append(image)
        cnt = 0
        x_signal_images = []
        for k in range(num_round):
            image = []
            for j in range(0,2*distance-1):
                row = []
                for i in range(0, 2*distance-1):
                    if (j,i) in mx_qubits:
                        row.append(signal_by_zerr[cnt])
                        cnt += 1
                    else:
                        row.append(0)
                image.append(row)
            x_signal_images.append(image)
        z_boundary_image = []
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                if j == 0 or j == 2*distance-2:
                    if i % 2 == 0:
                        row.append(1)
                    else:
                        row.append(0)
                else:
                    row.append(0)
            z_boundary_image.append(row)
        x_boundary_image = []
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                if i == 0 or i == 2*distance-2:
                    if j % 2 == 0:
                        row.append(1)
                    else:
                        row.append(0)
                else:
                    row.append(0)
            x_boundary_image.append(row)
        x_err_image = []
        cnt = 0
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                if (j,i) in data_qubits:
                    row.append(x_noise[cnt])
                    cnt += 1
                else:
                    row.append(0)
            x_err_image.append(row)
        z_err_image = []
        cnt = 0
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                if (j,i) in data_qubits:
                    row.append(z_noise[cnt])
                    cnt += 1
                else:
                    row.append(0)
            z_err_image.append(row)
        z_mea_image = []
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                row.append(0)
            z_mea_image.append(row)
        x_mea_image = []
        for j in range(0,2*distance-1):
            row = []
            for i in range(0, 2*distance-1):
                row.append(0)
            x_mea_image.append(row)
        samples.append(z_signal_images+x_signal_images+[z_boundary_image, x_boundary_image])
        labels.append([x_err_image, z_err_image, z_mea_image, x_mea_image])
    return np.array(samples), np.array(labels)