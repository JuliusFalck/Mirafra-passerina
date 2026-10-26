# pip install torch torchaudio librosa soundfile numpy scipy
import torch, torchaudio, librosa, numpy as np
from scipy.spatial.distance import cdist

TARGET_SR = 32000
WINDOW_SEC = 5.0

def load_audio_32k(path, target_sr=TARGET_SR):
    y, sr = librosa.load(path, sr=None, mono=True)
    if sr != target_sr:
        y = librosa.resample(y, orig_sr=sr, target_sr=target_sr)
    return y

def pad_or_trim(y, target_len):
    if len(y) >= target_len:
        return y[:target_len]
    z = np.zeros(target_len, dtype=y.dtype)
    z[:len(y)] = y
    return z

def embed_with_perch2(wave_32k: np.ndarray, sr=TARGET_SR) -> np.ndarray:
    """
    Replace this stub with the actual Perch 2 inference.
    Depending on the checkpoint, either:
      - feed raw audio and let the model do the spectrogram, or
      - compute log-mel here and feed that.
    Expected to return a 1D embedding vector (e.g., 1536-D).
    """
    # Example placeholder:
    # model = load_perch2_model_once()  # cache globally
    # with torch.no_grad():
    #     tensor = torch.tensor(wave_32k).float().unsqueeze(0)  # [1, T]
    #     emb = model.forward_audio(tensor, sr)                 # [1, D]
    # return emb.squeeze(0).cpu().numpy()
    raise NotImplementedError("Wire this to your Perch 2 checkpoint")

def song_embedding(path):
    y = load_audio_32k(path)
    target_len = int(WINDOW_SEC * TARGET_SR)
    y = pad_or_trim(y, target_len)
    return embed_with_perch2(y)

def pairwise_cosine(embeddings: np.ndarray):
    # cosine distance -> convert to similarity
    D = cdist(embeddings, embeddings, metric="cosine")
    S = 1.0 - D
    return S  # in [-1, 1]

# Usage:
# files = ["a.wav", "b.wav", "c.wav"]
# E = np.vstack([song_embedding(f) for f in files])  # [N, D]
# S = pairwise_cosine(E)  # similarity matrix
