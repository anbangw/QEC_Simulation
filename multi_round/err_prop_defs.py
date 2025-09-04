from enum import Enum

class ERR:
    def __init__(self, a):
        self.val = a
    def __eq__(self, other):
        return self.val == other.val
    def __mul__(self, other):
        if self.val == other.val:
            return ERR('I')
        elif self.val == 'I':
            return ERR(other.val)
        elif other.val == 'I':
            return ERR(self.val)
        else:
            for i in ['X', 'Y', 'Z']:
                if i != self.val and i != other.val:
                    return ERR(i)
    @staticmethod
    def I():
        return ERR('I')
    @staticmethod
    def X():
        return ERR('X')
    @staticmethod
    def Z():
        return ERR('Z')
    @staticmethod
    def Y():
        return ERR('Y')
    def __repr__(self):
        return self.val

def initial_no_error(qlist):
    init_err = {}
    for qb in qlist:
        init_err[qb] = ERR.I()
    return init_err

###### HELPER FUNCTIONS FOR SAMPLING NOISY CIRCUITS ######
# Sampling noisy Steane code stabilizer measurement circuits.
# recording signals after whole error detection circuit execution
import numpy as np
def propagate_through_noisy_reset(target_qubits, err_prob, init_err): # reset to zero
    current_err = init_err
    for qb in target_qubits:
        current_err[qb].val = 'I'
    reset_error_channel = ['I', 'X']
    err_reset = np.random.choice(reset_error_channel, len(target_qubits), p=[1-2/3*err_prob]+[2/3*err_prob])
    for i, qb in enumerate(target_qubits):
        current_err[qb] *= ERR(err_reset[i])

def propagate_through_noisy_measurement(target_qubits, err_prob, init_err): # z basis measurement
    current_err = init_err
    mea_error_channel = ['I', 'X']
    err_mea = np.random.choice(mea_error_channel, len(target_qubits), p=[1-2/3*err_prob]+[2/3*err_prob])
    for i, qb in enumerate(target_qubits):
        current_err[qb] *= ERR(err_mea[i])

def propagate_through_noisy_i_gate(target_qubits, err_prob, init_err):
    current_err = init_err
    clifford_gate_1q_error_channel = ['I', 'X', 'Z', 'Y']
    err_clifford_1q = np.random.choice(clifford_gate_1q_error_channel, len(target_qubits), p=[1-err_prob]+[1/3*err_prob]*3)
    for i, qb in enumerate(target_qubits):
        current_err[qb] *= ERR(err_clifford_1q[i])

def propagate_through_noisy_h_gate(target_qubits, err_prob, init_err):
    current_err = init_err
    for qb in target_qubits: # error channel changed
        if current_err[qb].val == 'Z':
            current_err[qb].val = 'X'
        elif current_err[qb].val == 'X':
            current_err[qb].val = 'Z'
    clifford_gate_1q_error_channel = ['I', 'X', 'Z', 'Y']
    err_clifford_1q = np.random.choice(clifford_gate_1q_error_channel, len(target_qubits), p=[1-err_prob]+[1/3*err_prob]*3)
    for i, qb in enumerate(target_qubits):
        current_err[qb] *= ERR(err_clifford_1q[i])

def propagate_through_noisy_x_gate(target_qubits, err_prob, init_err):
    propagate_through_noisy_i_gate(target_qubits, err_prob, init_err)

def propagate_through_noisy_z_gate(target_qubits, err_prob, init_err):
    propagate_through_noisy_i_gate(target_qubits, err_prob, init_err)

def propagate_through_noisy_cx_gate(control_qubits, target_qubits, err_prob, init_err): # tested okay
    current_err = init_err
    propagated_error_by_cx = initial_no_error(current_err.keys())
    for i, qb in enumerate(control_qubits):
        if current_err[qb].val in ['X', 'Y']:
            propagated_error_by_cx[target_qubits[i]] *= ERR.X()
    for i, qb in enumerate(target_qubits):
        if current_err[qb].val in ['Z', 'Y']:
            propagated_error_by_cx[control_qubits[i]] *= ERR.Z()
    for qb in current_err.keys():
        current_err[qb] *= propagated_error_by_cx[qb]
    # Errors from applying CX gates
    clifford_gate_1q_error_channel = ['I', 'X', 'Z', 'Y']
    clifford_gate_2q_error_channel = [i+j for i in clifford_gate_1q_error_channel for j in clifford_gate_1q_error_channel]
    err_clifford_2q = np.random.choice(clifford_gate_2q_error_channel, len(control_qubits), p=[1-15/16*err_prob]+[1/16*err_prob]*15)
    # print(err_clifford_2q)
    for i in range(len(control_qubits)):
        current_err[control_qubits[i]] *= ERR(err_clifford_2q[i][0])
        current_err[target_qubits[i]] *= ERR(err_clifford_2q[i][1])

def circuit_error_propagation(circ_str, err_prob, init_err):
    current_err = init_err
    for s in circ_str.split('\n'):
        s = s.strip()
        if s == '':
            continue
        else:
            # print(s)
            params = s.split()
            if params[0] == 'R': # reset
                propagate_through_noisy_reset([int(i) for i in params[1:]], err_prob, current_err)
            elif params[0] == 'RI': # ideal reset
                propagate_through_noisy_reset([int(i) for i in params[1:]], 0, current_err)
            elif params[0] == 'M': # measurement
                propagate_through_noisy_measurement([int(i) for i in params[1:]], err_prob, current_err)
            elif params[0] == 'MI': # ideal measurement
                propagate_through_noisy_measurement([int(i) for i in params[1:]], 0, current_err)
            elif params[0] == 'I':
                propagate_through_noisy_i_gate([int(i) for i in params[1:]], err_prob, current_err)
            elif params[0] == 'H':
                propagate_through_noisy_h_gate([int(i) for i in params[1:]], err_prob, current_err)
            elif params[0] == 'Z':
                propagate_through_noisy_z_gate([int(i) for i in params[1:]], err_prob, current_err)
            elif params[0] == 'X':
                propagate_through_noisy_x_gate([int(i) for i in params[1:]], err_prob, current_err)
            elif params[0] == 'CX':
                propagate_through_noisy_cx_gate([int(i) for i in params[1::2]], [int(i) for i in params[2::2]], err_prob, current_err)
            elif params[0] == '#':
                continue
            else:
                print('Circuit error propagation. Unsupported command:', s)
        
'''
def before_stabilizer_circuit_deploarizing(target_qubits, err_prob, init_err):
    propagate_through_i_gate(target_qubits, err_prob, init_err)
'''

###### HELPER FUNCTION FOR ERROR REPLAY ######
class errorLocation:
    def __init__(self, annotation, eid=None, error=None): # annotation is a string
        self.annotation = annotation
        self.eid = eid # the circuit position the error happens
        if error != None:
            self.error = error
        else:
            if self.annotation[0:2] == 'CX':
                self.error = 'II' # record the current error state
            else:
                self.error = 'I'
    def copy(self):
        return errorLocation(self.annotation, eid=self.eid, error=self.error)
    def __repr__(self):
        return self.annotation
    def error_channels(self): 
        # abides the circuit deploarizing model. You can change it by overloadding this function
        if self.annotation[0:2] == 'CX':
            clifford_gate_1q_error_channel = ['I', 'X', 'Z', 'Y']
            return [i+j for i in clifford_gate_1q_error_channel for j in clifford_gate_1q_error_channel]
        elif self.annotation[0] in ['I', 'H']:
            return ['I', 'X', 'Z', 'Y']
        elif self.annotation[0] in ['R', 'M']:
            return ['I', 'X']
        else:
            return []
    def width(self):
        return len(self.error)
    def error_sample(self, err_prob, num_shots=1):
        if self.annotation[0:2] == 'CX':
            clifford_gate_1q_error_channel = ['I', 'X', 'Z', 'Y']
            clifford_gate_2q_error_channel = [i+j for i in clifford_gate_1q_error_channel \
                                              for j in clifford_gate_1q_error_channel]
            err_clifford_2q = np.random.choice(clifford_gate_2q_error_channel, num_shots, \
                                               p=[1-15/16*err_prob]+[1/16*err_prob]*15)
            return err_clifford_2q
        elif self.annotation[0] in ['I', 'H']:
            clifford_gate_1q_error_channel = ['I', 'X', 'Z', 'Y']
            err_clifford_1q = np.random.choice(clifford_gate_1q_error_channel, num_shots, \
                                               p=[1-err_prob]+[1/3*err_prob]*3)
            return err_clifford_2q
        elif self.annotation[0] == 'R':
            reset_error_channel = ['I', 'X']
            err_reset = np.random.choice(reset_error_channel, num_shots, p=[1-2/3*err_prob]+[2/3*err_prob])
            return err_reset
        elif self.annotation[0] == 'M':
            mea_error_channel = ['I', 'X']
            err_mea = np.random.choice(mea_error_channel, num_shots, p=[1-2/3*err_prob]+[2/3*err_prob])
            return err_mea
class errorMechanism:
    def __init__(self):
        self._locations = []
        self.eid = 0
    def add_error_location(self, err_loc):
        self._locations.append(err_loc)
    def current_location(self):
        return self._locations[self.eid]
    def step_forward(self, step=1):
        self.eid = self.eid + step
    def step_reset(self):
        self.eid = 0
    def reset(self):
        self.eid = 0
        for loc in self._locations:
            loc.error = len(loc.error)*'I'
    def select_error_location(self, annotation, nth=1):
        for loc in self._locations:
            if loc.annotation == annotation:
                nth -= 1
                if nth == 0:
                    return loc
    def find_error_location(self, annotation, nth=1):
        for i, loc in enumerate(self._locations):
            if loc.annotation == annotation:
                nth -= 1
                if nth == 0:
                    return i
    def __len__(self):
        return len(self._locations)
    def __setitem__(self, eid, err_loc):
          self._locations[eid] = err_loc
    def __getitem__(self, eid):
          return self._locations[eid]
def propagate_through_reset_replay(target_qubits, init_err, error_mechanism): # reset to zero
    current_err = init_err
    for qb in target_qubits:
        current_err[qb].val = 'I'
    for qb in target_qubits:
        current_err[qb].val = error_mechanism.current_location().error
        error_mechanism.step_forward()
def propagate_through_measurement_replay(target_qubits, init_err, error_mechanism): # z basis measurement
    current_err = init_err
    for qb in target_qubits:
        current_err[qb] *= ERR(error_mechanism.current_location().error)
        error_mechanism.step_forward()
def propagate_through_i_gate_replay(target_qubits, init_err, error_mechanism):
    current_err = init_err
    for qb in target_qubits:
        current_err[qb] *= ERR(error_mechanism.current_location().error)
        error_mechanism.step_forward()
def propagate_through_x_gate_replay(target_qubits, init_err, error_mechanism):
    propagate_through_i_gate_replay(target_qubits, init_err, error_mechanism)
def propagate_through_z_gate_replay(target_qubits, init_err, error_mechanism):
    propagate_through_i_gate_replay(target_qubits, init_err, error_mechanism)
def propagate_through_h_gate_replay(target_qubits, init_err, error_mechanism):
    current_err = init_err
    for qb in target_qubits: # error channel changed
        if current_err[qb].val == 'Z':
            current_err[qb].val = 'X'
        elif current_err[qb].val == 'X':
            current_err[qb].val = 'Z'
    for qb in target_qubits:
        current_err[qb] *= ERR(error_mechanism.current_location().error)
        error_mechanism.step_forward()
def propagate_through_cx_gate_replay(control_qubits, target_qubits, init_err, error_mechanism):
    current_err = init_err
    propagated_error_by_cx = initial_no_error(current_err.keys())
    for i, qb in enumerate(control_qubits):
        if current_err[qb].val in ['X', 'Y']:
            propagated_error_by_cx[target_qubits[i]] *= ERR.X()
    for i, qb in enumerate(target_qubits):
        if current_err[qb].val in ['Z', 'Y']:
            propagated_error_by_cx[control_qubits[i]] *= ERR.Z()
    for qb in current_err.keys():
        current_err[qb] *= propagated_error_by_cx[qb]
    for i in range(len(control_qubits)):
        e1, e2 = error_mechanism.current_location().error
        current_err[control_qubits[i]] *= ERR(e1)
        current_err[target_qubits[i]] *= ERR(e2)
        error_mechanism.step_forward()
def circuit_error_propagation_replay(circ_str, init_err, error_mechanism):
    current_err = init_err
    for s in circ_str.split('\n'):
        s = s.strip()
        if s == '':
            continue
        else:
            params = s.split()
            if params[0] == 'R': # reset
                propagate_through_reset_replay([int(i) for i in params[1:]], current_err, error_mechanism)
            elif params[0] == 'RI':
                propagate_through_reset_replay([int(i) for i in params[1:]], current_err, error_mechanism)
            elif params[0] == 'M':
                propagate_through_measurement_replay([int(i) for i in params[1:]], current_err, error_mechanism)
            elif params[0] == 'MI':
                propagate_through_measurement_replay([int(i) for i in params[1:]], current_err, error_mechanism)
            elif params[0] == 'I':
                propagate_through_i_gate_replay([int(i) for i in params[1:]], current_err, error_mechanism)
            elif params[0] == 'H':
                propagate_through_h_gate_replay([int(i) for i in params[1:]], current_err, error_mechanism)
            elif params[0] == 'Z':
                propagate_through_z_gate_replay([int(i) for i in params[1:]], current_err, error_mechanism)
            elif params[0] == 'X':
                propagate_through_x_gate_replay([int(i) for i in params[1:]], current_err, error_mechanism)
            elif params[0] == 'CX':
                propagate_through_cx_gate_replay([int(i) for i in params[1::2]], [int(i) for i in params[2::2]], \
                                                 current_err, error_mechanism)
            elif params[0] == '#':
                continue
            else:
                print('Circuit error propagation replay. Unsupported command:', s)
def get_error_mechanism(circ_str):
    error_mechanism = errorMechanism()
    for s in circ_str.split('\n'):
        s = s.strip()
        if s == '':
            continue
        else:
            params = s.split()
            if params[0] in ['R', 'RI', 'M', 'MI', 'I', 'H', 'X', 'Z', 'T', 'S']:
                for qb in params[1:]:
                    error_mechanism.add_error_location(errorLocation(params[0]+' '+qb))
            elif params[0] == 'CX':
                ctr = params[1::2]
                tgt = params[2::2]
                for i in range(len(ctr)):
                    error_mechanism.add_error_location(errorLocation(params[0]+' '+ctr[i]+' '+tgt[i]))
            elif params[0] == '#':
                continue
            else:
                print('Extract error mechanism. Unsupported command:', s)
    return error_mechanism

###### SAMPLING FUNCTIONS ######
def list_kron(l1, l2): return [[i, j] for i in l1 for j in l2]

'''
Give error location, return error resulted from target X error locations.
'''
def sample_x_error_and_z_syndrome(code_circ, data_qubits, mx_qubits, mz_qubits, rounds=3,\
                                  err_location_samples=[], error_mechanism=None):
    if error_mechanism == None:
        error_mechanism = get_error_mechanism(code_circ*rounds)
    x_noise_list = []
    z_syndrome_list = []
    qubit_error_track = []
    loc_error_track = []
    x_error_channel = [['X'], ['XI', 'IX', 'XX']]
    # cnt = 0
    for loc in err_location_samples: # loc is a list of ErrorLocation
        if loc == []:
            error_mechanism.reset()
            init_err = initial_no_error(data_qubits+mx_qubits+mz_qubits)
            results = []
            err_res = []
            for i in range(rounds):
                circuit_error_propagation_replay(code_circ, init_err, error_mechanism)
                z_mea = [1 if init_err[qb].val in ['X', 'Y'] else 0 for qb in mz_qubits]
                x_noise = [1 if init_err[qb].val in ['X', 'Y'] else 0 for qb in data_qubits]
                results.append([z_mea, x_noise])
                err_res.append([init_err[e].val for e in init_err])
            x_noise_list.append(results[0][1])
            z_syndrome_list.append([signal[0] for signal in results])
            qubit_error_track.append(err_res)
            loc_error_track.append([])
            # cnt += 1
        else:
            err_list = []
            for l in loc: # l is a ErrorLocation object
                err_list.append(x_error_channel[error_mechanism[l].width()-1])
            all_loc_errors = [[]]
            for i in range(len(loc)):
                res = list_kron(all_loc_errors, [[err] for err in err_list[i]])
                all_loc_errors = [ec[0]+ec[1] for ec in res]
            # cnt += len(all_loc_errors)
            for error_chain in all_loc_errors:
                assert(len(error_chain) == len(loc))
                init_err = initial_no_error(data_qubits+mx_qubits+mz_qubits)
                results = []
                err_res = []
                loc_res = []
                error_mechanism.reset()
                for i, err in enumerate(error_chain):
                    error_mechanism[loc[i]].error = err
                    loc_res.append(error_mechanism[loc[i]].copy())
                for i in range(rounds):
                    circuit_error_propagation_replay(code_circ, init_err, error_mechanism)
                    z_mea = [1 if init_err[qb].val in ['X', 'Y'] else 0 for qb in mz_qubits]
                    x_noise = [1 if init_err[qb].val in ['X', 'Y'] else 0 for qb in data_qubits]
                    results.append([z_mea, x_noise])
                    err_res.append([init_err[e].val for e in init_err])
                x_noise_list.append(results[0][1])
                z_syndrome_list.append([signal[0] for signal in results])
                qubit_error_track.append(err_res)
                loc_error_track.append(loc_res)
    return z_syndrome_list, x_noise_list, qubit_error_track, loc_error_track
'''
Give error location, return error resulted from target Z error locations.
'''
def sample_z_error_and_x_syndrome(code_circ, data_qubits, mx_qubits, mz_qubits, rounds=3,\
                                  err_location_samples=[], error_mechanism=None):
    if error_mechanism == None:
        error_mechanism = get_error_mechanism(code_circ*rounds)
    z_noise_list = []
    x_syndrome_list = []
    qubit_error_track = []
    loc_error_track = []
    z_error_channel = [['Z'], ['ZI', 'IZ', 'ZZ']]
    # cnt = 0
    for loc in err_location_samples: # loc is a list of ErrorLocation
        if loc == []:
            error_mechanism.reset()
            init_err = initial_no_error(data_qubits+mx_qubits+mz_qubits)
            results = []
            err_res = []
            for i in range(rounds):
                circuit_error_propagation_replay(code_circ, init_err, error_mechanism)
                x_mea = [1 if init_err[qb].val in ['X', 'Y'] else 0 for qb in mx_qubits]
                z_noise = [1 if init_err[qb].val in ['Z', 'Y'] else 0 for qb in data_qubits]
                results.append([x_mea, z_noise])
                err_res.append([init_err[e].val for e in init_err])
            z_noise_list.append(results[0][1])
            x_syndrome_list.append([signal[0] for signal in results])
            qubit_error_track.append(err_res)
            loc_error_track.append([])
            # cnt += 1
        else:
            err_list = []
            for l in loc: # l is a ErrorLocation object
                err_list.append(z_error_channel[error_mechanism[l].width()-1])
            all_loc_errors = [[]]
            for i in range(len(loc)):
                res = list_kron(all_loc_errors, [[err] for err in err_list[i]])
                all_loc_errors = [ec[0]+ec[1] for ec in res]
            # cnt += len(all_loc_errors)
            for error_chain in all_loc_errors:
                assert(len(error_chain) == len(loc))
                init_err = initial_no_error(data_qubits+mx_qubits+mz_qubits)
                results = []
                err_res = []
                loc_res = []
                error_mechanism.reset()
                for i, err in enumerate(error_chain):
                    error_mechanism[loc[i]].error = err
                    loc_res.append(error_mechanism[loc[i]].copy())
                for i in range(rounds):
                    circuit_error_propagation_replay(code_circ, init_err, error_mechanism)
                    x_mea = [1 if init_err[qb].val in ['X', 'Y'] else 0 for qb in mx_qubits]
                    z_noise = [1 if init_err[qb].val in ['Z', 'Y'] else 0 for qb in data_qubits]
                    results.append([x_mea, z_noise])
                    err_res.append([init_err[e].val for e in init_err])
                z_noise_list.append(results[0][1])
                x_syndrome_list.append([signal[0] for signal in results])
                qubit_error_track.append(err_res)
                loc_error_track.append(loc_res)
    return x_syndrome_list, z_noise_list, qubit_error_track, loc_error_track
'''
Continuous generating error samples along the noisy circuit
'''
def continuous_circuit_error_random_sample(code_circ, data_qubits, mx_qubits, mz_qubits, err_prob, rounds=3, init_err=None):
    # data_qubits, mx_qubits, mz_qubits are the numbering of associated qubits
    if init_err == None:
        init_err = initial_no_error(data_qubits+mx_qubits+mz_qubits)
    results = []
    for i in range(rounds):
        circuit_error_propagation(code_circ, err_prob, init_err)
        x_mea_signal = [1 if init_err[qb].val in ['X', 'Y'] else 0 for qb in mx_qubits]
        z_mea_signal = [1 if init_err[qb].val in ['X', 'Y'] else 0 for qb in mz_qubits]
        x_noise = [1 if init_err[qb].val in ['X', 'Y'] else 0 for qb in data_qubits]
        z_noise = [1 if init_err[qb].val in ['Z', 'Y'] else 0 for qb in data_qubits]
        results.append([x_mea_signal, z_noise, z_mea_signal, x_noise])
    return results

'''
Given Z syndrome, check associated x error decoding ambiguity
'''
def check_error_location_detectable(syndrome_list, noise_list, observable):
    equ_err_loc = []
    for i in range(len(syndrome_list)):
        flag = False
        for equ_class in equ_err_loc:
            if syndrome_list[i] == syndrome_list[equ_class[0]]:
                equ_class.append(i)
                flag = True
        if flag == False:
            equ_err_loc.append([i])
    cnt = 0
    for equ_class in equ_err_loc:
        for i in equ_class:
            for j in equ_class:
                # Stabilizer commute not checked; only logical Z commute is checked
                corrected = (np.array(noise_list[i]) + np.array(noise_list[j])) % 2
                if sum(corrected*observable) % 2 != 0:
                    cnt += 1
    return cnt == 0 # True means no ambiguity