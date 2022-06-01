import functools
import multiprocessing
import os.path
import numpy as np
import pandas as pd
from scipy import signal, stats
from tqdm import tqdm
from more_itertools import consecutive_groups
from winrna_lib.helpers import Helpers

np.seterr(divide="ignore")


class WindowPeaks:
    def __init__(
        self,
        raw_signal: np.array,
        min_peak_distance,
        min_height,
        min_step_factor,
        is_reversed,
        prefix="",
    ):
        self.raw_signal = raw_signal
        self.windows = self.get_slicing_indexes(self.raw_signal)
        self.min_peak_distance = int(min_peak_distance)
        self.min_height = float(min_height)
        self.min_step_factor = float(min_step_factor)
        self.is_reversed = is_reversed
        self.prefix = prefix
        self.peaks_arr = self.call_signal_peaks()

    def call_signal_peaks(self) -> np.array:
        sig_deriv = (
            np.flipud(np.diff(np.flipud(self.raw_signal)))
            if self.is_reversed
            else np.diff(self.raw_signal)
        )
        all_peaks = []

        for window in tqdm(self.windows, desc="==> Calling peaks: ",
                           postfix="", bar_format='{desc} |{bar:20}| {percentage:3.0f}%'):
            window_sig_slice = sig_deriv[window[0]: window[1]]
            window_peaks = self._call_peaks_in_slice(
                window_sig_slice, self.min_peak_distance
            )

            if window_peaks is None:
                continue
            window_peaks[:, 0] += window[0]
            if window_peaks.size > 0:
                all_peaks.append(window_peaks)
        if not all_peaks:  # Exit if there is no peaks predicted
            return None
        all_peaks = np.concatenate(all_peaks, axis=0)
        if self.is_reversed:
            # all_peaks[:, 0] -= 1
            pass
        else:
            all_peaks[:, 0] += 1
        raw_heights = [self.raw_signal[x] for x in all_peaks[:, 0].astype(int)]
        mean_step_before = np.array([np.mean(self.raw_signal[x - 3: x]) for x in all_peaks[:, 0].astype(int)])
        mean_step_after = np.array([np.mean(self.raw_signal[x + 1: x + 4]) for x in all_peaks[:, 0].astype(int)])
        #step_factor = (mean_step_before / mean_step_after if self.is_reversed else mean_step_after / mean_step_before)
        step_factor = (mean_step_before / mean_step_after if self.is_reversed else mean_step_after / mean_step_before)
        # fold_change = np.abs(np.log2(mean_step_after / mean_step_before))\
        #    if self.is_reversed else \
        #    np.abs(np.log2(mean_step_before / mean_step_after))

        # plateau_height calculated as the average coverage for 30nt of peak height ignored if it contains zeros
        plateau_width = 12  # the optimum value is the minimum mapping length
        plateau_cov = np.array([self.raw_signal[x - plateau_width -1 : x] if self.is_reversed
                                else self.raw_signal[x : x + plateau_width]
                                for x in all_peaks[:, 0].astype(int)])

        mean_plateau_height = np.array([np.round(np.mean(pc), 2) for pc in plateau_cov])
        # np.round(fold_change, 2),
        all_peaks = np.stack(
            (
                all_peaks[:, 0],
                np.round(all_peaks[:, 1], 2),
                np.array(np.round(raw_heights, 2)),
                np.round(mean_step_before, 2),
                np.round(mean_step_after, 2),
                np.round(step_factor, 2),
                np.array(mean_plateau_height),
            ),
            axis=-1,
        )
        plateau_cov_filter__mask = np.array([np.all(pc) for pc in plateau_cov]) # filter interrupted signals
        raw_heights_filter_mask = all_peaks[:, 2] >= self.min_height
        step_factor_filter_mask = all_peaks[:, 5] >= self.min_step_factor
        all_peaks = all_peaks[plateau_cov_filter__mask & raw_heights_filter_mask & step_factor_filter_mask]  # apply filters
        return all_peaks

    @staticmethod
    def get_slicing_indexes(full_signal: np.array, slice_by=0.0, min_len=30):
        full_signal = np.stack(
            (np.array(list(range(0, full_signal.shape[0]))), full_signal), axis=-1
        )
        sliced_signal = full_signal[full_signal[:, 1] != slice_by]
        slice_indexes = [
            list(group) for group in consecutive_groups(sliced_signal[:, 0])
        ]
        return [
            (int(min(group)), int(max(group)))
            for group in slice_indexes
            if len(group) >= min_len
        ]

    @staticmethod
    def variable_iqr_threshold(sig_deriv):
        points = np.abs(sig_deriv[sig_deriv != 0])
        if points.size == 0:
            return None
        perc_list = [WindowPeaks._calc_custom_iqr(points, i) for i in range(75, 101, 1)]
        return perc_list[np.argmax(np.diff(perc_list)) - 1]

    @staticmethod
    def _calc_custom_iqr(data, prc):
        return stats.iqr(data, rng=(100 - prc, prc)) * 1.5 + np.percentile(data, prc)

    """
    @staticmethod
    def get_threshold_by_recursive_iqr(train_set, factor, sig_len, win_len=75):
        threshold_func = lambda data, factor: (np.percentile(data, 75) + stats.iqr(data) * 1.5) * factor
        thres = threshold_func(train_set, factor)
        if train_set[train_set >= thres].size > sig_len / win_len:
            factor += 0.1
            return WindowPeaks.get_threshold_by_recursive_iqr(train_set, factor, sig_len, win_len)
        else:
            return thres

    @staticmethod
    def get_threshold_by_zscore(train_set, min_score: float) -> np.array or None:
        sorted_train_set = np.sort(train_set)
        diff_train_set = np.diff(sorted_train_set)
        zscores = stats.zscore(diff_train_set)
        filt_mask = np.argwhere(zscores >= min_score)
        significant_values = sorted_train_set[filt_mask]
        if significant_values.shape[0] > 0:
            return np.min(significant_values)
    """
    @staticmethod
    def _call_peaks_in_slice(signal_slice: np.array, min_peak_distance=10):
        """
        # Core calling method
        threshold_func = lambda data: (np.percentile(data, 75) + stats.iqr(data)) * 1.5
        # Call peaks
        threshold_train_set = np.abs(signal_slice[signal_slice != 0])
        if threshold_train_set.size == 0:
            return None
        # threshold = WindowPeaks.get_threshold_by_recursive_iqr(threshold_train_set, 1, signal_slice.size)
        # threshold = threshold_func(threshold_train_set)
        threshold = WindowPeaks.get_threshold_by_zscore(threshold_train_set, 3)
        """
        threshold = WindowPeaks.variable_iqr_threshold(signal_slice)
        if threshold in [0, None]:
            return None
        peaks, peaks_props = signal.find_peaks(
            signal_slice, height=(threshold, None), distance=min_peak_distance)
        if peaks.size == 0:
            return None
        return np.stack((peaks, peaks_props["peak_heights"]), axis=-1)

    def get_peaks_df(self):
        # f"{self.prefix}_background_fold_change"
        if self.peaks_arr is None:
            return None
        peaks_df = pd.DataFrame(
            data=self.peaks_arr,
            columns=[
                "peak_index",
                f"{self.prefix}_diff_height",
                f"{self.prefix}_height",
                f"{self.prefix}_upstream",
                f"{self.prefix}_downstream",
                f"{self.prefix}_step_factor",
                f"{self.prefix}_mean_plateau_height",
            ],
            index=list(range(self.peaks_arr.shape[0])),
        )
        peaks_df["peak_index"] = peaks_df["peak_index"].astype(int)
        return peaks_df

    def get_bed_str(self, seqid: str):
        gff_df = self.get_peaks_df()
        if gff_df is None:
            return None
        columns = gff_df.columns.tolist()
        gff_df["seqid"] = seqid
        gff_df["start"] = gff_df["peak_index"] + 1
        gff_df["end"] = gff_df["peak_index"] + 1
        gff_df.drop(["peak_index"], inplace=True, axis=1)
        columns.remove("peak_index")
        gff_df["attributes"] = ""
        for indx in gff_df.index:
            gff_df.at[indx, "attributes"] = f"{self.prefix}_id={self.prefix}{indx}"
        for column in columns:
            gff_df["attributes"] += ";" + column + "=" + gff_df[column].astype(str)
        gff_df.drop(columns, inplace=True, axis=1)
        return gff_df.to_csv(index=False, sep="\t", header=False)

    def export_to_gff(
        self, out_path: str, seqid: str, strand: str, anno_source="_", anno_type="_"
    ):
        gff_columns = [
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ]
        gff_df = self.get_peaks_df()
        non_gff_columns = [x for x in gff_df.columns.tolist() if x not in gff_columns]
        gff_df["peak_index"] += 1
        gff_df.rename({"peak_index": "start"}, inplace=True, axis=1)
        non_gff_columns.remove("peak_index")
        gff_df["end"] = gff_df["start"]
        gff_df["seqid"] = seqid
        gff_df = Helpers.get_gff_df(
            gff_df,
            anno_source=anno_source,
            anno_type=anno_type,
            strand=strand,
            new_id=True,
        )
        gff_df.to_csv(os.path.abspath(out_path), index=False, sep="\t", header=False)
        print("GFF exported")
