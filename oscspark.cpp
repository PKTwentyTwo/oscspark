// Program for serching for spark + region -> region reactions.
#include "lifelib/pattern2.h"
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdint>
#include <set>
// Change the below variables:
std::string rule = "b3s23-a5";
std::string spark = "2o!";
std::string region = "3o$o$2o!";
std::string symmetry = "C1";
// Don't change these:
apg::lifetree<uint32_t, 1> lt(1000);
apg::lifetree<uint32_t, 3> lh(500);
std::vector<std::string> orientations = {"none", "rot90", "rot180", "rot270", "flip_x", "flip_y", "swap_xy", "swap_xy_flip"};
std::set<uint64_t> octodigests;
// Needed to get an RLE as a string:
std::string getrle(apg::pattern pat){
    std::ostringstream ss;
    pat.write_rle(ss);
    std::string s = ss.str();
    return s;
}
// Transform a pattern relative to the origin, depending on the symmetry.
apg::pattern transformpt(apg::pattern pt) {
    if (symmetry == "C2_1") {
        return pt + pt.transform("rot180", 0, 0);
    }
    if (symmetry == "C2_2") {
        return pt + pt.transform("rot180", 1, 0);
    }
    if (symmetry == "C2_4") {
        return pt + pt.transform("rot180", 1, 1);
    }
    return pt;
}
uint64_t GetPatternDigest(apg::pattern pat) {
    if (pat.empty()) { return 0; }
    int64_t bbox[4];
    pat.getrect(bbox);
    return pat.shift(0 - bbox[0], 0 - bbox[1]).digest();
}

// lifelib really should have an octodigest function at the C++ level.
// I should write a pull request later.
uint64_t GetPatternOctodigest(apg::pattern pat) {
    if (pat.empty()) { return 0; }

    uint64_t digests[8];
    uint8_t perms[8] = {27, 177, 78, 216, 39, 114, 141, 228};

    for (int i = 0; i < 8; i++) {

        auto lab = pat.getlab();

        apg::pattern tp(lab, lab->transform_recurse(pat.gethnode(), perms[i]), pat.getrule());

        int64_t bbox[4];
        tp.getrect(bbox);
        digests[i] = tp.shift(0 - bbox[0], 0 - bbox[1]).digest();
    }

    for (int i = 7; i > 0; i--) {
        for (int j = 0; j < i; j++) {
            if (digests[j] > digests[j+1]) {
                uint64_t c = digests[j];
                digests[j] = digests[j+1];
                digests[j+1] = c;
            }
        }
    }

    uint64_t acc = digests[0];
    uint64_t idx = 1;

    for (int i = 1; i < 8; i++) {
        if (digests[i] != digests[i-1]) {
            idx += 2;
            acc += digests[i] * idx;
        }
    }

    return acc;
}
// Borrowed from lifelib.cpp, and modified to remove the pointer nonsense.
apg::pattern* FindConnectedComponent(apg::pattern p_seed, apg::pattern p_backdrop, apg::pattern p_corona) {

    apg::pattern agglom = (p_seed) & (p_backdrop);
    apg::pattern cluster = agglom;

    while (cluster.nonempty()) {
        cluster = cluster.convolve(p_corona);
        cluster &= (p_backdrop);
        cluster -= agglom;
        agglom += cluster;
    }

    return new(std::nothrow) apg::pattern(agglom);
}
// Checks if a pattern regenerates:
int checkregen(apg::pattern pt, apg::pattern target, int gens, int ptgen) {
    int32_t i;
    apg::pattern ptcopy = pt.shift(0, 0);
    uint64_t targetdigest = GetPatternOctodigest(target);
    //std::cout << targetdigest << std::endl;
    for (i = 0; i < gens; i++) {
        uint64_t ptdigest = GetPatternOctodigest(ptcopy);
        // std::cout << targetdigest << " " << ptdigest << std::endl;
        if (ptdigest == targetdigest) {
            //If the spark(s) overlap with the region(s), this will catch them by returning 0.
            return i + ptgen;
        }
        if (ptcopy.empty()) {
            break;
        }
        ptcopy = ptcopy[1];
    }
    //No match; return 0 to signify a failure:
    return 0;
}
//Function to do the actual searching:
void searchpt(int gens) {
    apg::pattern sparkpt(&lt, spark, rule);
    apg::pattern regionpt(&lt, region, rule);
    std::vector<apg::pattern> sparkpts;
    sparkpts.reserve(8);
    int i, j;
    for (auto i : orientations) {
        apg::pattern tpt = sparkpt.shift(0, 0);
        if (i != "none") {
            tpt = sparkpt.transform(i, 0, 0);
        }
        uint64_t digest = GetPatternDigest(tpt);
        if (octodigests.find(digest) == octodigests.end()) {
            octodigests.insert(digest);
            sparkpts.push_back(sparkpt.transform(i, 0, 0));
        }
    }
    apg::pattern halo(&lt, "3o$3o$3o!", rule);
    halo = halo.centre();
    //Stamp collection is 
    apg::pattern stamp(&lh, "b!", rule+"History");
    std::cout << "Searching " << sparkpts.size() << " orientations of spark..." << std::endl;
    int patterns = 0;
    for (auto sparkpt : sparkpts) {    
        int gen = 0;
        int ptgen = 0;
        while (gen < gens) {
            apg::pattern ptcopy1 = regionpt[gen];
            int64_t bbox[4] = {0, 0, 0, 0};
            ptcopy1.getrect(bbox);
            for (i = bbox[0]-6; i < bbox[0] + bbox[2]+5; i++) {
                for (j = bbox[1]-6; j < bbox[1] + bbox[3]+5; j++) {
                    //std::cout << i << "," << j << std::endl;
                    apg::pattern ptcopy2 = ptcopy1[0];
                    ptcopy2 += sparkpt.shift(i, j);
                    
                    ptgen = checkregen(ptcopy2, regionpt, 150, gen);
                    if (ptgen > 2) {
                        uint64_t octodigest = GetPatternOctodigest(ptcopy2);
                        if (octodigests.find(octodigest) == octodigests.end()) { 
                            if (FindConnectedComponent(sparkpt.shift(i, j).onecell(), ptcopy2, halo)->totalPopulation() == sparkpt.totalPopulation()) {
                                //std::cout << "Regenerates after " << ptgen << " generations:" << std::endl;
                                //std::cout << getrle(ptcopy2) << std::endl;
                                stamp += ptcopy2.shift(100 * (patterns % 10), 100 * (patterns / 10));
                                octodigests.insert(octodigest);
                                patterns++;
                                if (patterns%5 == 0) {
                                    std::cout << patterns << " partials found." << std::endl;
                                }
                            }
                        }
                    }
                }
            }
            gen += 1;
            if (gen%10 == 0) {
                std::cout << gen << " generations searched." << std::endl;
            }
            if (ptcopy1.empty()) {
                std::cout << "Region dies out at generation " << gen << std::endl;
                break;
            }
        }
    }
    std::cout << getrle(stamp) << std::endl;
}
        
int main(int argc, char* argv[]) {
    if (argc < 4) {
        if (argc == 2) {
            if (strcmp(argv[1], "-h") == 0)  {
                std::cout << "Usage: oscspark <rule> 'region_rle' 'spark_rle'\nrule:     The rule being searched. Run './compile.sh <rule>' first.\n\nregion_rle:   The RLE of the active region you want to spark.\nMust be headerless and enclosed by quotes.\n\nspark_rle: The RLE of the spark you want to use.\nMust be headerless and enclosed by quotes - e.g 'o!'" << std::endl;
                return 0;
            }
        }
        std::cout << "Usage: oscspark <rule> 'region_rle' 'spark_rle'\nTry './oscspark -h' for more information." << std::endl;
        return 0;
    }
    rule = argv[1];
    region = argv[2];
    spark = argv[3];
    searchpt(100);
}
