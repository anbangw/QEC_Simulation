UP = lambda x: (x[0]-1, x[1])
DW = lambda x: (x[0]+1, x[1])
LF = lambda x: (x[0], x[1]-1)
RF = lambda x: (x[0], x[1]+1)
def list_xor(a, b):
    return [a[i]^b[i] for i in range(len(a))]
import pymatching
from math import log
import numpy as np
from .err_prop_defs import *
class surfaceCode:
    def __init__(self, distance):
        self.data_qubits = []
        for j in range(0,2*distance-1):
            for i in range(0, 2*distance-1):
                if j % 2 == i % 2:
                    # j: row id, i: column id
                    self.data_qubits.append((j, i))
        self.mz_qubits = []
        for j in range(0,2*distance-1):
            for i in range(0, 2*distance-1):
                if i % 2 == 1 and j % 2 == 0:
                    self.mz_qubits.append((j, i))
        self.mx_qubits = []
        for j in range(0,2*distance-1):
            for i in range(0, 2*distance-1):
                if j % 2 == 1 and i % 2 == 0:
                    self.mx_qubits.append((j, i))
        self.distance = distance
        self.qubit_ordering = {}
        cnt = 0
        for qb in self.data_qubits+self.mz_qubits+self.mx_qubits:
            self.qubit_ordering[qb] = cnt
            cnt += 1
    def get_z_parity_matrix(self):
        parity_matrix = []
        for zq in self.mz_qubits:
            if zq[0] == 0: # upper boundary
                row = []
                for qb in self.data_qubits:
                    if qb in [LF(zq), RF(zq), DW(zq)]:
                        row.append(1)
                    else:
                        row.append(0)
                parity_matrix.append(row)
            elif zq[0] == 2*self.distance-1-1:
                row = []
                for qb in self.data_qubits:
                    if qb in [LF(zq), RF(zq), UP(zq)]:
                        row.append(1)
                    else:
                        row.append(0)
                parity_matrix.append(row)
            else:
                row = []
                for qb in self.data_qubits:
                    if qb in [LF(zq), RF(zq), UP(zq), DW(zq)]:
                        row.append(1)
                    else:
                        row.append(0)
                parity_matrix.append(row)
        return parity_matrix
    def get_x_parity_matrix(self):
        parity_matrix = []
        for xq in self.mx_qubits:
            if xq[1] == 0: # upper boundary
                row = []
                for qb in self.data_qubits:
                    if qb in [RF(xq), UP(xq), DW(xq)]:
                        row.append(1)
                    else:
                        row.append(0)
                parity_matrix.append(row)
            elif xq[1] == 2*self.distance-1-1:
                row = []
                for qb in self.data_qubits:
                    if qb in [LF(xq), UP(xq), DW(xq)]:
                        row.append(1)
                    else:
                        row.append(0)
                parity_matrix.append(row)
            else:
                row = []
                for qb in self.data_qubits:
                    if qb in [LF(xq), RF(xq), UP(xq), DW(xq)]:
                        row.append(1)
                    else:
                        row.append(0)
                parity_matrix.append(row)
        return parity_matrix
    def get_x_observable(self):
        observable = []
        for qb in self.data_qubits:
            if qb[0] == self.distance - 1:
                observable.append(1)
            else:
                observable.append(0)
        return observable
    def get_z_observable(self):
        observable = []
        for qb in self.data_qubits:
            if qb[1] == self.distance - 1:
                observable.append(1)
            else:
                observable.append(0)
        return observable
    # detect X error using pymatching
    def get_z_signal_matching_graph(self, err_prob=0.01, rounds=3, weight_func=None):
        mz_num = {}
        for i, mzq in enumerate(self.mz_qubits):
            mz_num[mzq] = i
        dq_num = {}
        for i, dq in enumerate(self.data_qubits):
            dq_num[dq] = i
        if weight_func == None:
            weight_func = lambda p: log((1-p)/p)
        m = pymatching.Matching()
        for L in range(rounds):
            for mzq in self.mz_qubits:
                if mzq[1] == 1:
                    m.add_boundary_edge(mz_num[mzq]+L*len(self.mz_qubits), fault_ids=dq_num[LF(mzq)], \
                                        weight=weight_func(8/15*err_prob), merge_strategy='smallest-weight')
                elif mzq[1] == 2*self.distance-2-1:
                    m.add_boundary_edge(mz_num[mzq]+L*len(self.mz_qubits), fault_ids=dq_num[RF(mzq)], 
                                        weight=weight_func(8/15*err_prob), merge_strategy='smallest-weight')
                    # m.add_boundary_edge(mz_num[mzq], fault_ids=1, weight=1, error_probability=0.1)
            for mzq in self.mz_qubits:
                if mzq[1] != 2*self.distance-2-1:
                    m.add_edge(mz_num[mzq]+L*len(self.mz_qubits), mz_num[RF(RF(mzq))]+L*len(self.mz_qubits), \
                        fault_ids=dq_num[RF(mzq)], weight=weight_func(2*8/15*err_prob+4*4/15*err_prob+2/3*err_prob), \
                               merge_strategy='smallest-weight')
                    #, error_probability=0.05)
            for mzq in self.mz_qubits:
                if mzq[0] != 2*self.distance-2:
                    m.add_edge(mz_num[mzq]+L*len(self.mz_qubits), mz_num[DW(DW(mzq))]+L*len(self.mz_qubits), \
                        fault_ids=dq_num[DW(mzq)], weight=weight_func(2*4/15*err_prob+2/3*err_prob), 
                               merge_strategy='smallest-weight')
        for L in range(rounds-1):
            for mzq in self.mz_qubits:
                # no correction can be applied for measurement errors
                m.add_edge(mz_num[mzq]+L*len(self.mz_qubits), mz_num[mzq]+(L+1)*len(self.mz_qubits), \
                        fault_ids=None, weight=weight_func(4*4/15*err_prob+2/3*err_prob), 
                           merge_strategy='smallest-weight')
            for mzq in self.mz_qubits:
                if mzq[0] != 2*self.distance-2:
                    m.add_edge(mz_num[mzq]+L*len(self.mz_qubits), mz_num[DW(DW(mzq))]+(L+1)*len(self.mz_qubits), \
                        fault_ids=dq_num[DW(mzq)], weight=weight_func(4*4/15*err_prob), 
                               merge_strategy='smallest-weight')
            for mzq in self.mz_qubits:
                if mzq[1] != 2*self.distance-2-1:
                    m.add_edge(mz_num[mzq]+L*len(self.mz_qubits), mz_num[RF(RF(mzq))]+(L+1)*len(self.mz_qubits), \
                        fault_ids=dq_num[RF(mzq)], weight=weight_func(2*4/15*err_prob), 
                               merge_strategy='smallest-weight')
            for mzq in self.mz_qubits:
                if mzq[1] != 1 and mzq[0] != 2*self.distance-2:
                    m.add_edge(mz_num[mzq]+L*len(self.mz_qubits), mz_num[DW(DW(LF(LF(mzq))))]+(L+1)*len(self.mz_qubits), \
                        fault_ids={dq_num[DW(mzq)],dq_num[LF(DW(DW(mzq)))]}, weight=weight_func(2*4/15*err_prob), 
                               merge_strategy='smallest-weight')
        return m
    # detect z error using pymatching
    def get_x_signal_matching_graph(self, err_prob=0.01, rounds=3, weight_func=None):
        mx_num = {}
        for i, mxq in enumerate(self.mx_qubits):
            mx_num[mxq] = i
        dq_num = {}
        for i, dq in enumerate(self.data_qubits):
            dq_num[dq] = i
        if weight_func == None:
            weight_func = lambda p: log((1-p)/p)
        m = pymatching.Matching()
        for L in range(rounds):
            for mxq in self.mx_qubits:
                if mxq[0] == 1:
                    m.add_boundary_edge(mx_num[mxq]+L*len(self.mx_qubits), fault_ids=dq_num[UP(mxq)], \
                                        weight=weight_func(8/15*err_prob), merge_strategy='smallest-weight')
                elif mxq[0] == 2*self.distance-2-1:
                    m.add_boundary_edge(mx_num[mxq]+L*len(self.mx_qubits), fault_ids=dq_num[DW(mxq)], 
                                        weight=weight_func(8/15*err_prob), merge_strategy='smallest-weight')
            for mxq in self.mx_qubits:
                if mxq[1] != 2*self.distance-2:
                    m.add_edge(mx_num[mxq]+L*len(self.mx_qubits), mx_num[RF(RF(mxq))]+L*len(self.mx_qubits), \
                        fault_ids=dq_num[RF(mxq)], weight=weight_func(2*8/15*err_prob+4*4/15*err_prob+2/3*err_prob), \
                               merge_strategy='smallest-weight')
            for mxq in self.mx_qubits:
                if mxq[0] != 2*self.distance-2-1:
                    m.add_edge(mx_num[mxq]+L*len(self.mx_qubits), mx_num[DW(DW(mxq))]+L*len(self.mx_qubits), \
                        fault_ids=dq_num[DW(mxq)], weight=weight_func(2*4/15*err_prob+2/3*err_prob), 
                               merge_strategy='smallest-weight')
        for L in range(rounds-1):
            for mxq in self.mx_qubits:
                # no correction can be applied for measurement errors
                m.add_edge(mx_num[mxq]+L*len(self.mx_qubits), mx_num[mxq]+(L+1)*len(self.mx_qubits), \
                        fault_ids=None, weight=weight_func(4*4/15*err_prob+2/3*err_prob), 
                           merge_strategy='smallest-weight')
            for mxq in self.mx_qubits:
                if mxq[0] != 2*self.distance-2-1:
                    m.add_edge(mx_num[mxq]+L*len(self.mx_qubits), mx_num[DW(DW(mxq))]+(L+1)*len(self.mx_qubits), \
                        fault_ids=dq_num[DW(mxq)], weight=weight_func(4*4/15*err_prob), 
                               merge_strategy='smallest-weight')
            for mxq in self.mx_qubits:
                if mxq[1] != 2*self.distance-2:
                    m.add_edge(mx_num[mxq]+L*len(self.mx_qubits), mx_num[RF(RF(mxq))]+(L+1)*len(self.mx_qubits), \
                        fault_ids=dq_num[RF(mxq)], weight=weight_func(2*4/15*err_prob), 
                               merge_strategy='smallest-weight')
            for mxq in self.mx_qubits:
                if mxq[1] != 0 and mxq[0] != 2*self.distance-2-1:
                    m.add_edge(mx_num[mxq]+L*len(self.mx_qubits), mx_num[DW(DW(LF(LF(mxq))))]+(L+1)*len(self.mx_qubits), \
                        fault_ids={dq_num[DW(mxq)],dq_num[LF(DW(DW(mxq)))]}, weight=weight_func(2*4/15*err_prob), 
                               merge_strategy='smallest-weight')
        return m
    # get circuit string for sampling
    def __get_sampling_circ(self, rounds=1):
        circ = ""
        q_num = self.qubit_ordering
        circ += 'I'
        for qb in self.data_qubits:
            circ += f' {q_num[qb]}'
        circ += '\n'
        for i in range(rounds-1):
            circ += 'R'
            for qb in self.mz_qubits+self.mx_qubits:
                circ += f' {q_num[qb]}'
            circ += '\n'
            circ += 'I'
            for qb in self.mz_qubits:
                circ += f' {q_num[qb]}'
            circ += '\n'
            circ += 'H'
            for qb in self.mx_qubits:
                circ += f' {q_num[qb]}'
            circ += '\n'
            # step 1
            circ += 'CX'
            for qb in self.mz_qubits:
                if qb[0] != 0:
                    circ += f' {q_num[UP(qb)]} {q_num[qb]}'
            for qb in self.mx_qubits:
                circ += f' {q_num[qb]} {q_num[UP(qb)]}'
            circ += '\n'
            # step 2
            circ += 'CX'
            for qb in self.mz_qubits:
                circ += f' {q_num[LF(qb)]} {q_num[qb]}'
            for qb in self.mx_qubits:
                if qb[1] != 0:
                    circ += f' {q_num[qb]} {q_num[LF(qb)]}'
            circ += '\n'
            # step 3
            circ += 'CX'
            for qb in self.mz_qubits:
                circ += f' {q_num[RF(qb)]} {q_num[qb]}'
            for qb in self.mx_qubits:
                if qb[1] != 2*self.distance-2:
                    circ += f' {q_num[qb]} {q_num[RF(qb)]}'
            circ += '\n'
            # step 4
            circ += 'CX'
            for qb in self.mz_qubits:
                if qb[0] != 2*self.distance-2:
                    circ += f' {q_num[DW(qb)]} {q_num[qb]}'
            for qb in self.mx_qubits:
                circ += f' {q_num[qb]} {q_num[DW(qb)]}'
            circ += '\n'
            circ += 'I'
            for qb in self.mz_qubits:
                circ += f' {q_num[qb]}'
            circ += '\n'
            circ += 'H'
            for qb in self.mx_qubits:
                circ += f' {q_num[qb]}'
            circ += '\n'
            circ += 'M'
            for qb in self.mz_qubits+self.mx_qubits:
                circ += f' {q_num[qb]}'
            circ += '\n'
            circ += 'ROUND END\n'
        # final ideal measurement
        circ += 'IR'
        for qb in self.mz_qubits+self.mx_qubits:
            circ += f' {q_num[qb]}'
        circ += '\n'
        circ += 'II'
        for qb in self.mz_qubits:
            circ += f' {q_num[qb]}'
        circ += '\n'
        circ += 'IH'
        for qb in self.mx_qubits:
            circ += f' {q_num[qb]}'
        circ += '\n'
        # step 1
        circ += 'ICX'
        for qb in self.mz_qubits:
            if qb[0] != 0:
                circ += f' {q_num[UP(qb)]} {q_num[qb]}'
        for qb in self.mx_qubits:
            circ += f' {q_num[qb]} {q_num[UP(qb)]}'
        circ += '\n'
        # step 2
        circ += 'ICX'
        for qb in self.mz_qubits:
            circ += f' {q_num[LF(qb)]} {q_num[qb]}'
        for qb in self.mx_qubits:
            if qb[1] != 0:
                circ += f' {q_num[qb]} {q_num[LF(qb)]}'
        circ += '\n'
        # step 3
        circ += 'ICX'
        for qb in self.mz_qubits:
            circ += f' {q_num[RF(qb)]} {q_num[qb]}'
        for qb in self.mx_qubits:
            if qb[1] != 2*self.distance-2:
                circ += f' {q_num[qb]} {q_num[RF(qb)]}'
        circ += '\n'
        # step 4
        circ += 'ICX'
        for qb in self.mz_qubits:
            if qb[0] != 2*self.distance-2:
                circ += f' {q_num[DW(qb)]} {q_num[qb]}'
        for qb in self.mx_qubits:
            circ += f' {q_num[qb]} {q_num[DW(qb)]}'
        circ += '\n'
        circ += 'II'
        for qb in self.mz_qubits:
            circ += f' {q_num[qb]}'
        circ += '\n'
        circ += 'IH'
        for qb in self.mx_qubits:
            circ += f' {q_num[qb]}'
        circ += '\n'
        circ += 'IM'
        for qb in self.mz_qubits+self.mx_qubits:
            circ += f' {q_num[qb]}'
        circ += '\n'
        circ += 'ROUND END\n'
        return q_num, circ
    def _sample(self, rounds, err_prob, num_shots):
        q_num, circ = self.__get_sampling_circ(rounds)
        instructions = []
        for s in circ.split('\n'):
            s = s.strip()
            if s == '':
                continue
            else:
                # print(s)
                params = s.split()
                if params[0] in ['R', 'IR', 'M', 'IM', 'I', 'II', 'H', 'IH', 'Z', 'IZ', 'X', 'IX']: # reset
                    instructions.append([params[0], [int(i) for i in params[1:]]])
                elif params[0] == 'CX':
                    instructions.append([params[0], [int(i) for i in params[1::2]], [int(i) for i in params[2::2]]])
                elif params[0] == 'ICX':
                    instructions.append([params[0], [int(i) for i in params[1::2]], [int(i) for i in params[2::2]]])
                elif params[0] == 'ROUND':
                    instructions.append([params[0]])
                else:
                    continue
        # print(instructions)
        res = []
        for i in range(num_shots):
            mea_data = []
            init_err = {}
            for qb in self.data_qubits+self.mz_qubits+self.mx_qubits:
                init_err[q_num[qb]] = ERR.I()
            for ins in instructions:
                if ins[0] == 'R':
                    propagate_through_noisy_reset(ins[1], err_prob, init_err)
                elif ins[0] == 'M':
                    propagate_through_noisy_measurement(ins[1], err_prob, init_err)
                elif ins[0] == 'I':
                    propagate_through_noisy_i_gate(ins[1], err_prob, init_err)
                elif ins[0] == 'H':
                    propagate_through_noisy_h_gate(ins[1], err_prob, init_err)
                elif ins[0] == 'Z':
                    propagate_through_noisy_z_gate(ins[1], err_prob, init_err)
                elif ins[0] == 'X':
                    propagate_through_noisy_x_gate(ins[1], err_prob, init_err)
                elif ins[0] == 'CX':
                    propagate_through_noisy_cx_gate(ins[1], ins[2], err_prob, init_err)
                elif ins[0] == 'IR':
                    propagate_through_noisy_reset(ins[1], 0, init_err)
                elif ins[0] == 'IM':
                    propagate_through_noisy_measurement(ins[1], 0, init_err)
                elif ins[0] == 'II':
                    propagate_through_noisy_i_gate(ins[1], 0, init_err)
                elif ins[0] == 'IH':
                    propagate_through_noisy_h_gate(ins[1], 0, init_err)
                elif ins[0] == 'IZ':
                    propagate_through_noisy_z_gate(ins[1], 0, init_err)
                elif ins[0] == 'IX':
                    propagate_through_noisy_x_gate(ins[1], 0, init_err)
                elif ins[0] == 'ICX':
                    propagate_through_noisy_cx_gate(ins[1], ins[2], 0, init_err)
                else: # round ends
                    x_mea_signal = [1 if init_err[q_num[qb]].val in ['X', 'Y'] else 0 for qb in self.mx_qubits]
                    z_mea_signal = [1 if init_err[q_num[qb]].val in ['X', 'Y'] else 0 for qb in self.mz_qubits]
                    x_noise = [1 if init_err[q_num[qb]].val in ['X', 'Y'] else 0 for qb in self.data_qubits]
                    z_noise = [1 if init_err[q_num[qb]].val in ['Z', 'Y'] else 0 for qb in self.data_qubits]
                    mea_data.append([x_mea_signal, z_noise, z_mea_signal, x_noise])
            res.append(mea_data)
        return res
    # z error: 0, x error: 1
    def _sample_errors(self, rounds, err_prob, num_shots, err_type):
        test_data = []
        data = self._sample(rounds, err_prob, num_shots)
        syndrome_pos = 2*err_type
        for test_sample in data:
            detection_event = [k for k in test_sample[0][syndrome_pos]]
            for test_round in range(1, len(test_sample)):
                # detection event are triggered if any signal is flipped.
                detection_event += list_xor(test_sample[test_round-1][syndrome_pos], test_sample[test_round][syndrome_pos])
            noise = [k for k in test_sample[len(test_sample)-1][syndrome_pos+1]]
            test_data.append([detection_event, noise])
        return test_data
    def get_sampling_circ(self, rounds=1):
        return self.__get_sampling_circ(rounds)
    def sample_x_error(self, rounds, err_prob, num_shots):
        return self._sample_errors(rounds, err_prob, num_shots, 1)
    def sample_z_error(self, rounds, err_prob, num_shots):
        return self._sample_errors(rounds, err_prob, num_shots, 0)
    def x_error_decoding_failure(self, pred, noise):
        corrected = (np.array(pred) + np.array(noise)) % 2
        parity_matrix = np.array(self.get_z_parity_matrix())
        observable = np.array(self.get_z_observable())
        if np.allclose((parity_matrix@corrected)%2,0) and np.allclose((observable@corrected)%2,0):
            return 0
        else:
            return 1
    def z_error_decoding_failure(self, pred, noise):
        corrected = (np.array(pred) + np.array(noise)) % 2
        parity_matrix = np.array(self.get_x_parity_matrix())
        observable = np.array(self.get_x_observable())
        if np.allclose((parity_matrix@corrected)%2,0) and np.allclose((observable@corrected)%2,0):
            return 0
        else:
            return 1