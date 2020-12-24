# IMPORTANT, only tested using PYTHON 3!
import numpy as np
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import matplotlib.pylab as plt
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.special import logsumexp
import pandas as pd
import string

################
# note: if you are modifying the alphabet
# make sure last character is "-" (gap)
################
alphabet = "ARNDCQEGHILKMFPSTWYV-"
states = len(alphabet)

# map amino acids to integers (A->0, R->1, etc)
a2n = dict((a,n) for n,a in enumerate(alphabet))
aa2int = lambda x: a2n.get(x,a2n['-'])

# from fasta
def parse_fasta(filename, a3m=False):
  '''function to parse fasta file'''
  
  if a3m:
    # for a3m files the lowercase letters are removed
    # as these do not align to the query sequence
    rm_lc = str.maketrans(dict.fromkeys(string.ascii_lowercase))
    
  header, sequence = [],[]
  lines = open(filename, "r")
  for line in lines:
    line = line.rstrip()
    if line[0] == ">":
      header.append(line[1:])
      sequence.append([])
    else:
      if a3m: line = line.translate(rm_lc)
      else: line = line.upper()
      sequence[-1].append(line)
  lines.close()
  sequence = [''.join(seq) for seq in sequence]
  
  return np.array(header), np.array(sequence)
'''
def parse_fasta(filename):
  function to parse fasta
  header = []
  sequence = []
  lines = open(filename, "r")
  for line in lines:
    line = line.rstrip()
    if line[0] == ">":
      header.append(line[1:])
      sequence.append([])
    else:
      sequence[-1].append(line)
  lines.close()
  sequence = [''.join(seq) for seq in sequence]
  return np.array(header), np.array(sequence)
'''

def filt_gaps(msa, gap_cutoff=0.5):
  '''filters alignment to remove gappy positions'''
  frac_gaps = np.mean((msa == states-1).astype(np.float),0)
  non_gaps = np.where(frac_gaps < gap_cutoff)[0]
  return msa[:,non_gaps], non_gaps

def get_eff(msa, eff_cutoff=0.8):
  '''compute effective weight for each sequence'''  
  msa_sm = 1.0 - squareform(pdist(msa,"hamming"))
  msa_w = (msa_sm >= eff_cutoff).astype(np.float)
  msa_w = 1.0/np.sum(msa_w,-1)
  return msa_w

def str2int(x):
  '''convert a list of strings into list of integers'''
  # Example: ["ACD","EFG"] -> [[0,4,3], [6,13,7]]
  if x.dtype.type is np.str_:
    if x.ndim == 0: return np.array([aa2int(aa) for aa in x])
    else: return np.array([[aa2int(aa) for aa in seq] for seq in x])
  else:
    return x
  
def split_train_test(seqs, frac_test=0.1):
  # shuffle data
  x = np.copy(seqs)
  np.random.shuffle(x[1:])

  # fraction of data used for testing
  split = int(len(x) * (1.0-frac_test))

  # split training/test datasets
  return x[:split], x[split:]

def mk_msa(seqs, gap_cutoff=0.5, eff_cutoff=0.8):
  '''converts list of sequences to MSA (Multiple Sequence Alignment)'''
  # =============================================================================
  # The function takes a list of sequences (strings) and returns a (dict)ionary
  # containing the following:
  # =============================================================================
  # BEFORE GAP REMOVAL
  # -----------------------------------------------------------------------------
  # msa_ori   msa
  # ncol_ori  number of columns
  # -----------------------------------------------------------------------------
  # AFTER GAP REMOVAL
  # By default, columns with ≥ 50% gaps are removed. This makes things a
  # little complicated, as we need to keep track of which positions were removed.
  # -----------------------------------------------------------------------------
  # msa       msa
  # ncol      number of columns
  # v_idx     index of positions kept
  # -----------------------------------------------------------------------------
  # weights   weight for each sequence (based on sequence identity)
  # nrow      number of rows (sequences)
  # neff      number of effective sequences sum(weights)
  # =============================================================================
  
  msa_ori = str2int(seqs)

  # remove positions with more than > 50% gaps
  msa, v_idx = filt_gaps(msa_ori)
  
  # compute effective weight for each sequence
  msa_weights = get_eff(msa, eff_cutoff)
    
  return {"msa_ori":msa_ori,
          "msa":msa,
          "weights":msa_weights,
          "neff":np.sum(msa_weights),
          "v_idx":v_idx,
          "nrow":msa.shape[0],
          "ncol":msa.shape[1],
          "ncol_ori":msa_ori.shape[1]}

# optimizer
def opt_adam(loss, name, var_list=None, lr=1.0, b1=0.9, b2=0.999, b_fix=False):
  # adam optimizer
  # Note: this is a modified version of adam optimizer. More specifically, we replace "vt"
  # with sum(g*g) instead of (g*g). Furthmore, we find that disabling the bias correction
  # (b_fix=False) speeds up convergence for our case.
  
  if var_list is None: var_list = tf.trainable_variables() 
  gradients = tf.gradients(loss,var_list)
  if b_fix: t = tf.Variable(0.0,"t")
  opt = []
  for n,(x,g) in enumerate(zip(var_list,gradients)):
    if g is not None:
      ini = dict(initializer=tf.zeros_initializer,trainable=False)
      mt = tf.get_variable(name+"_mt_"+str(n),shape=list(x.shape), **ini)
      vt = tf.get_variable(name+"_vt_"+str(n),shape=[], **ini)
      
      mt_tmp = b1*mt+(1-b1)*g
      vt_tmp = b2*vt+(1-b2)*tf.reduce_sum(tf.square(g))
      lr_tmp = lr/(tf.sqrt(vt_tmp) + 1e-8)

      if b_fix: lr_tmp = lr_tmp * tf.sqrt(1-tf.pow(b2,t))/(1-tf.pow(b1,t))

      opt.append(x.assign_add(-lr_tmp * mt_tmp))
      opt.append(vt.assign(vt_tmp))
      opt.append(mt.assign(mt_tmp))
        
  if b_fix: opt.append(t.assign_add(1.0))
  return(tf.group(opt))

def GREMLIN(msa,
            opt_iter=100,
            opt_rate=1.0,
            batch_size=None,
            lam_v=0.01,
            lam_w=0.01,
            scale_lam_w=True,
            v=None,
            w=None,
            ignore_gap=True):
  
  '''fit params of MRF (Markov Random Field) given MSA (multiple sequence alignment)'''
  # ==========================================================================
  # this function takes a MSA (dict)ionary, from mk_msa() and returns a MRF
  # (dict)ionary containing the following:
  # ==========================================================================
  # len       full length
  # v_idx     index of positions (mapping back to full length)
  # v         2-body term
  # w         2-body term
  # ==========================================================================
  # WARNING: The mrf is over the msa after gap removal. "v_idx" and "len" are
  # important for mapping the MRF back to the original MSA.
  # ==========================================================================
  
  ########################################
  # SETUP COMPUTE GRAPH
  ########################################
  # reset tensorflow graph
  tf.reset_default_graph()
  
  # length of sequence
  ncol = msa["ncol"] 
  
  # input msa (multiple sequence alignment) 
  MSA = tf.placeholder(tf.int32,shape=(None,ncol),name="msa")
  
  # input msa weights
  MSA_weights = tf.placeholder(tf.float32, shape=(None,), name="msa_weights")
  
  # one-hot encode msa
  OH_MSA = tf.one_hot(MSA,states)
  
  if ignore_gap:
    ncat = states - 1
    NO_GAP = 1.0 - OH_MSA[...,-1] 
    OH_MSA = OH_MSA[...,:ncat]
    
  else:
    ncat = states
  
  ########################################
  # V: 1-body-term of the MRF
  ########################################
  V = tf.get_variable(name="V",
                          shape=[ncol,ncat],
                          initializer=tf.zeros_initializer)
  
  ########################################
  # W: 2-body-term of the MRF
  ########################################
  W_tmp = tf.get_variable(name="W",
                          shape=[ncol,ncat,ncol,ncat],
                          initializer=tf.zeros_initializer)  
  
  # symmetrize W
  W = W_tmp + tf.transpose(W_tmp,[2,3,0,1])
  
  # set diagonal to zero
  W = W * (1-np.eye(ncol))[:,None,:,None]

  ########################################
  # Pseudo-Log-Likelihood
  ########################################
  # V + W
  VW = V + tf.tensordot(OH_MSA,W,2)
  
  # hamiltonian
  H = tf.reduce_sum(OH_MSA*VW,-1)
  
  # local Z (parition function)
  Z = tf.reduce_logsumexp(VW,-1)

  PLL = H - Z
  if ignore_gap:
    PLL = PLL * NO_GAP  

  PLL = tf.reduce_sum(PLL,-1)  
  PLL = tf.reduce_sum(MSA_weights * PLL)/tf.reduce_sum(MSA_weights)

  ########################################
  # Regularization
  ########################################
  L2 = lambda x: tf.reduce_sum(tf.square(x))
  L2_V = lam_v * L2(V)
  L2_W = lam_w * L2(W) * 0.5
  
  if scale_lam_w:
    L2_W = L2_W * (ncol-1) * (states-1)
  
  ########################################
  # Loss Function
  ########################################
  # loss function to minimize
  loss = -PLL + (L2_V + L2_W) / msa["neff"]
  
  # optimizer
  opt = opt_adam(loss,"adam",lr=opt_rate)
  
  ########################################
  # Input Generator
  ########################################
  all_idx = np.arange(msa["nrow"])
  def feed(feed_all=False):
    if batch_size is None or feed_all:
      return {MSA:msa["msa"], MSA_weights:msa["weights"]}
    else:
      batch_idx = np.random.choice(all_idx,size=batch_size)
      return {MSA:msa["msa"][batch_idx], MSA_weights:msa["weights"][batch_idx]}
  
  ########################################
  # OPTIMIZE
  ########################################
  with tf.Session() as sess:
    
    # initialize variables V and W
    sess.run(tf.global_variables_initializer())

    # initialize V
    if v is None:
      oh_msa = np.eye(states)[msa["msa"]]
      if ignore_gap: oh_msa = oh_msa[...,:-1]
      
      pseudo_count = 0.01 * np.log(msa["neff"])
      f_v = np.einsum("nla,n->la",oh_msa,msa["weights"])
      V_ini = np.log(f_v + pseudo_count)
      if lam_v > 0:
        V_ini = V_ini - np.mean(V_ini,axis=-1,keepdims=True)
      sess.run(V.assign(V_ini))
      
    else:
      sess.run(V.assign(v))

    # initialize W
    if w is not None:
      sess.run(W_tmp.assign(w * 0.5))
      
    # compute loss across all data
    get_loss = lambda: np.round(sess.run(loss,feed(True)) * msa["neff"],2)

    print("starting",get_loss())      
    for i in range(opt_iter):
      sess.run(opt,feed())  
      if (i+1) % int(opt_iter/10) == 0:
        print("iter",(i+1),get_loss())
    
    # save the V and W parameters of the MRF
    V_ = sess.run(V)
    W_ = sess.run(W)
    
  ########################################
  # return MRF
  ########################################
  no_gap_states = states - 1
  mrf = {"v": V_[:,:no_gap_states],
         "w": W_[:,:no_gap_states,:,:no_gap_states],
         "v_idx": msa["v_idx"],
         "len": msa["ncol_ori"]}
  
  return mrf