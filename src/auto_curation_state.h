#ifndef AUTO_CURATION_STATE_H 
#define AUTO_CURATION_STATE_H

#include "utilsPretextView.h"
#include "showWindowData.h"
#include "genomeData.h"
#include <algorithm>
#include <vector>

namespace {

inline u32 median_u32_3(u32 a, u32 b, u32 c)
{
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    return b;
}

// Hi-C map contig ids can flip for single pixels at boundaries; median-of-3 along the 1D axis
// stabilizes runs so Pixel Sort does not build huge duplicate fragment lists (e.g. 130,131,130,...).
inline u32 smoothed_contig_id_1d(
    map_state* Map_State,
    u32 i,
    u32 start_pixel,
    u32 end_pixel)
{
    u32 left = (i > start_pixel) ? Map_State->contigIds[i - 1] : Map_State->contigIds[i];
    u32 mid = Map_State->contigIds[i];
    u32 right = (i < end_pixel) ? Map_State->contigIds[i + 1] : Map_State->contigIds[i];
    return median_u32_3(left, mid, right);
}

// Median-of-3 can still alternate every pixel when raw contigIds alternate. Collapse a trailing
// suffix that strictly alternates between exactly two ids (length >= 4) to [left, other] in map order.
inline void collapse_alternating_two_contig_suffix(std::vector<u32>& ids)
{
    const size_t n = ids.size();
    if (n < 4) return;

    const size_t e = n - 1;
    const u32 a = ids[e];
    const u32 b = ids[e - 1];
    if (a == b) return;

    size_t s = e - 1;
    while (s > 0)
    {
        const u32 x = ids[s - 1];
        if (x != a && x != b) break;
        if (x == ids[s]) break;
        --s;
    }

    const size_t len = n - s;
    if (len < 4) return;

    for (size_t i = s; i + 1 < n; ++i)
        if (ids[i] == ids[i + 1]) return;

    for (size_t i = s; i < n; ++i)
        if (ids[i] != a && ids[i] != b) return;

    const u32 left = ids[s];
    const u32 other = (left == a) ? b : a;
    ids.erase(ids.begin() + (ptrdiff_t)s, ids.end());
    ids.push_back(left);
    ids.push_back(other);
}

inline std::vector<u32> build_smoothed_selected_run_ids(
    map_state* Map_State,
    u32 sp,
    u32 ep)
{
    std::vector<u32> out;
    if (sp > ep) return out;
    u32 last = smoothed_contig_id_1d(Map_State, sp, sp, ep);
    out.push_back(last);
    for (u32 i = sp + 1; i <= ep; i++)
    {
        u32 cur = smoothed_contig_id_1d(Map_State, i, sp, ep);
        if (cur != last)
        {
            out.push_back(cur);
            last = cur;
        }
    }
    collapse_alternating_two_contig_suffix(out);
    return out;
}

} // namespace

struct SelectArea
{   
public:
    SelectArea() {};
    u08 select_flag = 0; // no selected area set
    u32 start_pixel = 0;
    u32 end_pixel = 0;
    s32 source_frag_id = -1;
    s32 sink_frag_id = -1;
    std::vector<s32> selected_frag_ids;


    void clearup()
    {
        this->select_flag = 0;
        this->start_pixel = 0;
        this->end_pixel = 0;
        this->source_frag_id = -1;
        this->sink_frag_id = -1;
        this->selected_frag_ids.clear();
    }

    /*
    len in pixel, including the source and sink fragments
    */
    u32 get_selected_len(
        const map_contigs* Contigs,
        const bool used_for_cluster_flag=false
    ) const
    {
        u32 len = 0;
        if (this->source_frag_id>=0 && !used_for_cluster_flag)
        {
            len += Contigs->contigs_arr[this->source_frag_id].length;
        }
        for (const auto& frag_id : this->selected_frag_ids)
        {
            len += Contigs->contigs_arr[frag_id].length;
        }
        if (this->sink_frag_id>=0 && !used_for_cluster_flag)
        {
            len += Contigs->contigs_arr[this->sink_frag_id].length;
        }
        return len;
    }

    /*
    number of fragments to be sorted, including the source and sink fragments
    */
    u32 get_to_sort_frags_num()
    {
        return (this->source_frag_id>=0) + this->selected_frag_ids.size() +  (this->sink_frag_id>=0); 
    }

    /*
    id of fragments to be sorted, including the source and sink fragments
    */
    std::vector<u32> get_to_sort_frags_id(const map_contigs* Contigs)
    {
        std::vector<u32> to_sort_frags;
        if (this->source_frag_id>=0 )
        {
            to_sort_frags.push_back(this->source_frag_id);
        }
        for (const auto& frag_id : this->selected_frag_ids)
        {
            to_sort_frags.push_back(frag_id);
        }
        if (this->sink_frag_id>=0 )
        {
            to_sort_frags.push_back(this->sink_frag_id);
        }
        return to_sort_frags;
    }
};


struct AutoCurationState
{

private:
    s32 start_pixel = -1; // select contigs for sort
    s32 end_pixel = -1;   // select contigs for sort

public:
    f32 mask_color[4] = {0.906f, 0.03921f, 0.949f, 0.2f}; // used to show the selected area.
    f32 mask_color_default[4] = {0.906f, 0.03921f, 0.949f, 0.2f}; // default color
    u08 selected_or_exclude = 0; // 0 for select, 1 for exclude
    u08 select_mode = 0; // 0: deactive, 1: set mode, 2: delete mode from start, 3: delete mode from end
    // Variables to store the current values
    u32 smallest_frag_size_in_pixel = 2;
    f32 link_score_threshold = 0.4f;

    f32 auto_cut_threshold = 0.05f;
    u32 auto_cut_diag_window_for_pixel_mean = 8;
    u32 auto_cut_smallest_frag_size_in_pixel = 8;

    // cluster according to the hap name
    u08 hap_cluster_flag = 0; 

    // kmeans cluster
    u32 num_clusters = 1;
    const s32 min_frag_num_for_cluster = 4;

    // Variables for the editing UI state
    bool show_autoSort_erase_confirm_popup = false;
    bool show_autoSort_redo_confirm_popup = false;
    //sorting mode
    u32 sort_mode = 1; // 0: union find, 1: fuse union find, 2 deep fuse, 3 yahs
    std::vector<std::string> sort_mode_names = {"Union Find", "Fuse", "Deep Fuse"};
    
    // pixel sort
    char frag_size_buf[16];
    char score_threshold_buf[16];

    // pixel cut
    char auto_cut_threshold_buf[16];
    char auto_cut_diag_window_for_pixel_mean_buf[16];
    char auto_cut_smallest_frag_size_in_pixel_buf[16];
    bool auto_cut_with_extension = true;  // considering the gap extension while cutting
    s32 auto_cut_gap_loc_threshold = 3;   // if the distance between the calculated cut loc and the gap is less than this value, the cut will be made at the gap

    AutoCurationState()
    {
        set_buf();
    }
    
    void set_buf()
    {   
        // sort: smallest_frag_size
        memset(this->frag_size_buf, 0, sizeof(frag_size_buf));
        snprintf((char*)frag_size_buf, sizeof(frag_size_buf), "%u", this->smallest_frag_size_in_pixel);
        // sort: link_score_threshold
        memset(score_threshold_buf, 0, sizeof(score_threshold_buf));
        snprintf((char*)this->score_threshold_buf, sizeof(this->score_threshold_buf), "%.3f", this->link_score_threshold);

        // // cut: auto_cut_threshold
        // memset(this->auto_cut_threshold_buf, 0, sizeof(this->auto_cut_threshold_buf));
        // snprintf((char*)this->auto_cut_threshold_buf, sizeof(this->auto_cut_threshold_buf), "%.3f", this->auto_cut_threshold);
        // // cut: auto_cut_diag_window_for_pixel_mean
        // memset(this->auto_cut_diag_window_for_pixel_mean_buf, 0, sizeof(this->auto_cut_diag_window_for_pixel_mean_buf));
        // snprintf((char*)this->auto_cut_diag_window_for_pixel_mean_buf, sizeof(this->auto_cut_diag_window_for_pixel_mean_buf), "%u", this->auto_cut_diag_window_for_pixel_mean);
        // // cut: auto_cut_smallest_frag_size_in_pixel
        // memset(this->auto_cut_smallest_frag_size_in_pixel_buf, 0, sizeof(this->auto_cut_smallest_frag_size_in_pixel_buf));
        // snprintf((char*)this->auto_cut_smallest_frag_size_in_pixel_buf, sizeof(this->auto_cut_smallest_frag_size_in_pixel_buf), "%u", this->auto_cut_smallest_frag_size_in_pixel);


        fmt::format_to((char*)auto_cut_threshold_buf, "{:.3f}", this->auto_cut_threshold);
        fmt::format_to((char*)auto_cut_diag_window_for_pixel_mean_buf, "{}", this->auto_cut_diag_window_for_pixel_mean);
        fmt::format_to((char*)auto_cut_smallest_frag_size_in_pixel_buf, "{}", this->auto_cut_smallest_frag_size_in_pixel);
    }

    void set_mask_color(f32* color)
    {
        for (u32 i = 0; i < 4; i++)
        {
            this->mask_color[i] = color[i];
        }
    }

    void set_default_mask_color()
    {
        for (u32 i = 0; i < 4; i++)
        {
            this->mask_color[i] = this->mask_color_default[i];
        }
    }

    void update_value_from_buf()
    {   
        // sort 
        this->smallest_frag_size_in_pixel = (u32)atoi((char*)this->frag_size_buf);
        this->link_score_threshold = (float)atof((char*)this->score_threshold_buf);
        if (this->link_score_threshold > 1.0f || this->link_score_threshold < 0.0f) 
        {   
            printf("[Auto Sort] Warning: link Score Threshold should be in the range of [0, 1]\n");
            this->link_score_threshold = std::max(0.0f, std::min(1.0f, this->link_score_threshold));
        }

        // cut
        this->auto_cut_threshold = (f32)atof((char*)this->auto_cut_threshold_buf);
        this->auto_cut_diag_window_for_pixel_mean = (u32)atoi((char*)this->auto_cut_diag_window_for_pixel_mean_buf);
        this->auto_cut_smallest_frag_size_in_pixel = (u32)atoi((char*)this->auto_cut_smallest_frag_size_in_pixel_buf);

        this->set_buf();
    }

    std::string get_sort_mode_name() const
    {
        return sort_mode_names[sort_mode];
    }

    void clear() 
    {
        this->start_pixel = -1;
        this->end_pixel = -1;
    }

    void set_start(u32 start_pix) 
    {
        this->start_pixel = start_pix;
        this->select_mode = 2;
    }

    void set_end(u32 end_pix)
    {
        this->end_pixel = end_pix;
        this->select_mode = 3;
    }

    s32 get_start()
    {   
        return this->start_pixel;
    }

    s32 get_end()
    {
        return this->end_pixel;
    }

    u32 get_selected_fragments_num(
        map_state* Map_State,
        u32 number_of_pixels_1D
    ) 
    {
        // Bounds: Map_State->contigIds is [0, number_of_pixels_1D - 1]. Reject >= number_of_pixels_1D
        // (a strict `>` check would wrongly allow end_pixel == number_of_pixels_1D and overrun the array).
        if (number_of_pixels_1D == 0 ||
            this->start_pixel >= (s32)number_of_pixels_1D ||
            this->end_pixel >= (s32)number_of_pixels_1D ||
            this->start_pixel < 0 ||
            this->end_pixel < 0)
        {
            return 0;
        }
        if (this->start_pixel > this->end_pixel)
        {
            std::swap(this->start_pixel, this->end_pixel);
        }

        u32 sp = (u32)this->start_pixel, ep = (u32)this->end_pixel;
        return (u32)build_smoothed_selected_run_ids(Map_State, sp, ep).size();
    }

    std::vector<u32> get_selected_fragments_id(
        map_state* Map_State,
        u32 number_of_pixels_1D
    )
    {
        std::vector<u32> selected_frag_ids={};
        // Same bounds as get_selected_fragments_num / get_selected_fragments.
        if (number_of_pixels_1D == 0 ||
            this->start_pixel >= (s32)number_of_pixels_1D ||
            this->end_pixel >= (s32)number_of_pixels_1D ||
            this->start_pixel < 0 ||
            this->end_pixel < 0)
        {
            return selected_frag_ids;
        }
        if (this->start_pixel > this->end_pixel)
        {
            std::swap(this->start_pixel, this->end_pixel);
        }

        u32 sp = (u32)this->start_pixel, ep = (u32)this->end_pixel;
        return build_smoothed_selected_run_ids(Map_State, sp, ep);
    }


    void get_selected_fragments(
        SelectArea& select_area, 
        map_state* Map_State, 
        u32 number_of_pixels_1D, 
        map_contigs* Contigs,
        u08 exclude_flag=false
    ) // return cluster_flag
    {   
        // Walk [start_pixel, end_pixel] and record one fragment id per contig run (when contigIds changes).
        // - Bounds: indices must satisfy 0 <= start/end < number_of_pixels_1D (use >= below; do not use
        //   end_pixel > number_of_pixels_1D only, or end_pixel == number_of_pixels_1D can index past the array).
        // - clearup() every time: avoids stale source_frag_id / sink_frag_id when the new selection has no
        //   neighbour on the left (start at 0) or right (end at last pixel); those stay -1 by design.
        if (number_of_pixels_1D == 0 ||
            this->start_pixel >= (s32)number_of_pixels_1D ||
            this->end_pixel >= (s32)number_of_pixels_1D ||
            this->start_pixel < 0 ||
            this->end_pixel < 0)
        {   
            return ;
        }
        if (this->start_pixel > this->end_pixel)
        {
            std::swap(this->start_pixel, this->end_pixel);
        }

        select_area.clearup();

        // Copy run list into select_area (use this->start_pixel / this->end_pixel consistently for the loop).
        select_area.start_pixel = (u32)this->start_pixel;
        select_area.end_pixel = (u32)this->end_pixel;
        u32 sp = select_area.start_pixel, ep = select_area.end_pixel;
        for (u32 cid : build_smoothed_selected_run_ids(Map_State, sp, ep))
            select_area.selected_frag_ids.push_back((s32)cid);

        // Source = fragment on the pixel immediately before the selection (contigIds[start_pixel - 1]).
        // If start_pixel == 0 there is no left neighbour; source_frag_id remains -1 (callers must handle).
        select_area.source_frag_id = -1;
        if (this->start_pixel > 0)
        {   
            u32 frag_id = Map_State->contigIds[(u32)this->start_pixel - 1];
            if (Contigs->contigs_arr[frag_id].length >= this->smallest_frag_size_in_pixel)
            {
                select_area.source_frag_id = (s32)frag_id;
            }
            else
            {   
                fmt::print("The source_frag_len ({}) < smallest_length ({}), not set the source frag.\n", Contigs->contigs_arr[frag_id].length, this->smallest_frag_size_in_pixel);
            }
        }
        
        // Sink = fragment on the pixel immediately after the selection (contigIds[end_pixel + 1]).
        // If end_pixel is the last map pixel, end_pixel + 1 is out of range; sink_frag_id remains -1.
        select_area.sink_frag_id = -1;
        if ((u32)this->end_pixel + 1u < number_of_pixels_1D)
        {   
            u32 frag_id = Map_State->contigIds[(u32)this->end_pixel + 1u];
            if (Contigs->contigs_arr[frag_id].length >= this->smallest_frag_size_in_pixel)
            {
                select_area.sink_frag_id = (s32)frag_id;
            }
            else
            {   
                fmt::print("The tail_frag_len ({}) < smallest_length ({}), not set the tail frag.\n", Contigs->contigs_arr[frag_id].length, this->smallest_frag_size_in_pixel);
            }
        }
        // set the select_area to valid
        if (!select_area.selected_frag_ids.empty()) select_area.select_flag = 1;

        return;
    }

    void update_sort_area(
        u32 x_pixel, map_state* Map_State, u32 num_pixels_1d
    )
    {   
        if (this->select_mode == 1) // expand area
        {   
            if (this->start_pixel == -1) this->start_pixel = x_pixel;
            if (this->end_pixel == -1) this->end_pixel = x_pixel;
            if (x_pixel < this->start_pixel)
            {
                this->start_pixel = x_pixel;
            }
            else if (x_pixel > this->end_pixel)
            {
                this->end_pixel = x_pixel;
            }
            else if (this->start_pixel < x_pixel && x_pixel < this->end_pixel)
            {   
                u32 dis_to_start = std::abs(this->start_pixel - (s32)x_pixel);
                u32 dis_to_end = std::abs(this->end_pixel - (s32)x_pixel);
                if (dis_to_start < dis_to_end)
                {
                    this->start_pixel = x_pixel;
                }
                else 
                {
                    this->end_pixel = x_pixel;
                }
            }
        }
        else if (this->select_mode == 2) // narrow from start
        {
            if (x_pixel > this->start_pixel && x_pixel <= this->end_pixel)
            {
                this->start_pixel = x_pixel;
            }
            else if (x_pixel> end_pixel)
            {
                clear();
            }
        }
        else if (this->select_mode == 3) // narrow from end
        {
            if (this->start_pixel <= x_pixel && x_pixel < this->end_pixel)
            {
                this->end_pixel = x_pixel;
            }
            else if (x_pixel < start_pixel)
            {
                clear();
            }
        }

        // update start_pixel and end_pixel to cover the range
        if (select_mode == 1 && this->start_pixel!=-1 && this->end_pixel!=-1)
        {
            u32 tmp_pixel = this->start_pixel;
            while ( this->start_pixel>0 && Map_State->contigIds[tmp_pixel] == Map_State->contigIds[this->start_pixel])
            {
                this->start_pixel -- ;
            }
            if (Map_State->contigIds[tmp_pixel] != Map_State->contigIds[this->start_pixel])  this->start_pixel ++ ;
            tmp_pixel = this->end_pixel;
            while ( (this->end_pixel < num_pixels_1d - 1) && Map_State->contigIds[tmp_pixel] == Map_State->contigIds[this->end_pixel])
            {
                this->end_pixel ++ ;
            }
            if (Map_State->contigIds[tmp_pixel] != Map_State->contigIds[this->end_pixel]) this->end_pixel --;
        }
        return;
    }

    void adjust_cut_threshold(s32 delta)
    {
        this->auto_cut_threshold += delta * 0.002f;
        this->auto_cut_threshold = std::max(0.0f, std::min(1.0f, this->auto_cut_threshold));
        fmt::format_to((char*)auto_cut_threshold_buf, "{:.4f}", this->auto_cut_threshold);
    }

    void change_sort_mode(s32 delta)
    {
        this->sort_mode = (this->sort_mode + this->sort_mode_names.size() + (delta>0?1:-1)) % this->sort_mode_names.size();
    }

};

#endif // AUTO_CURATION_STATE_H 