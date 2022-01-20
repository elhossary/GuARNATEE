import os.path
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy import signal, stats


class WindowPeaks:

    def __init__(self, raw_signal: np.array, windows: dict, min_peak_distance: int,
                 threshold_factor: float, min_height, is_reversed):
        self.raw_signal = raw_signal
        self.windows = windows
        self.min_peak_distance = min_peak_distance
        self.threshold_factor = threshold_factor
        self.min_height = min_height
        self.is_reversed = is_reversed
        self.peaks_arr = self.call_signal_peaks()

    def call_signal_peaks(self) -> np.array:
        sig_deriv = np.flipud(np.diff(np.flipud(self.raw_signal))) if self.is_reversed \
            else np.diff(self.raw_signal)
        all_peaks = []
        for group in tqdm(self.windows.keys(), desc="Progress: "):
            window = self.windows[group]
            window_sig_slice = sig_deriv[window[0]: window[1]]
            window_peaks = self._call_peaks_in_slice(window_sig_slice,
                                                     self.min_peak_distance,
                                                     self.threshold_factor)
            if window_peaks is None:
                continue
            window_peaks[:, 0] += window[0]
            all_peaks.append(window_peaks)
        all_peaks = np.concatenate(all_peaks, axis=0)
        if self.is_reversed:
            #all_peaks[:, 0] -= 1
            pass
        else:
            all_peaks[:, 0] += 1
        raw_heights = [self.raw_signal[x] for x in all_peaks[:, 0].astype(int)]
        all_peaks = np.stack((all_peaks[:, 0], all_peaks[:, 1], np.array(raw_heights)), axis=-1)
        all_peaks = all_peaks[all_peaks[:, 2] >= self.min_height]
        return all_peaks

    @staticmethod
    def _call_peaks_in_slice(signal_slice, min_peak_distance, threshold_factor):
        peaks, peaks_props = signal.find_peaks(signal_slice,
                                               width=(None, None),
                                               distance=min_peak_distance,
                                               height=(None, None),
                                               threshold=(None, None),
                                               prominence=(None, None),
                                               rel_height=0.5,
                                               plateau_size=(1, None))
        if peaks.size == 0:
            return None

        # Generate threshold
        q3 = np.percentile(peaks_props["peak_heights"], 75)
        iqr = stats.iqr(peaks_props["peak_heights"])
        threshold = q3 + (iqr * threshold_factor)
        # Filter by threshold
        peaks_arr = np.stack((peaks, peaks_props["peak_heights"]), axis=-1)
        return peaks_arr[peaks_arr[:, 1] >= threshold]

    def get_peaks_df(self):
        peaks_df = pd.DataFrame(data=self.peaks_arr,
                                columns=["peak_index", "diff_height", "height"],
                                index=list(range(self.peaks_arr.shape[0])))
        peaks_df["peak_index"] = peaks_df["peak_index"].astype(int)
        return peaks_df

    def get_bed_str(self, seqid: str):
        gff_df = self.get_peaks_df()
        columns = gff_df.columns.tolist()
        gff_df["seqid"] = seqid
        gff_df["start"] = gff_df["peak_index"] + 1
        gff_df["end"] = gff_df["peak_index"] + 1
        gff_df.drop(columns, inplace=True, axis=1)
        return gff_df.to_csv(index=False, sep="\t", header=False)

    def export_to_gff(self, out_path: str, seqid: str, strand: str, anno_source="_", anno_type="_"):
        gff_columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        gff_df = self.get_peaks_df()
        non_gff_columns = [x for x in gff_df.columns.tolist() if x not in gff_columns]
        gff_df["peak_index"] += 1
        gff_df.rename({"peak_index": "start"}, inplace=True, axis=1)
        non_gff_columns.remove("peak_index")
        gff_df["end"] = gff_df["start"]
        gff_df["seqid"] = seqid
        gff_df["source"] = anno_source
        gff_df["type"] = anno_type
        gff_df["score"] = "."
        gff_df["strand"] = strand
        gff_df["phase"] = "."
        for i in gff_df.index:
            gff_df.at[i, "attributes"] = f'ID={anno_type}_{gff_df.at[i, "seqid"]}{strand}_{i}'\
                                         f';name={anno_type}_{gff_df.at[i, "seqid"]}{strand}_{i}'
            for col in non_gff_columns:
                gff_df.at[i, "attributes"] += f';{col}={gff_df.at[i, col]}'

        gff_df.drop(non_gff_columns, inplace=True, axis=1)
        gff_df = gff_df.reindex(columns=gff_columns)
        gff_df.to_csv(os.path.abspath(out_path), index=False, sep="\t", header=False)
        print("GFF exported")