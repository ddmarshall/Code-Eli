/*********************************************************************************
* Copyright (c) 2014 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#ifndef piecewise_cst_airfoil_creator_test_suite_hpp
#define piecewise_cst_airfoil_creator_test_suite_hpp

#include <cmath>    // std::pow, std::exp

#include <typeinfo> // typeid
#include <string>   // std::string
#include <sstream>  // std::stringstream
#include <iomanip>  // std::setw
#include <limits>   // std::numeric_limits

#include "eli/constants/math.hpp"
#include "eli/mutil/fd/d1o2.hpp"
#include "eli/mutil/fd/d2o2.hpp"
#include "eli/geom/curve/piecewise.hpp"
#include "eli/geom/curve/piecewise_cst_airfoil_creator.hpp"
#include "eli/geom/curve/piecewise_cst_airfoil_fitter.hpp"
#include "eli/geom/curve/pseudo/cst_airfoil.hpp"

template<typename data__>
class piecewise_cst_airfoil_creator_test_suite : public Test::Suite
{
  private:
    typedef eli::geom::curve::piecewise<eli::geom::curve::bezier, data__, 3> piecewise_curve_type;
    typedef typename piecewise_curve_type::curve_type curve_type;
    typedef typename piecewise_curve_type::point_type point_type;
    typedef typename piecewise_curve_type::control_point_type control_point_type;
    typedef typename piecewise_curve_type::data_type data_type;
    typedef typename piecewise_curve_type::index_type index_type;
    typedef typename piecewise_curve_type::tolerance_type tolerance_type;
    typedef eli::geom::curve::pseudo::cst_airfoil<data_type> cst_airfoil_type;
    typedef typename cst_airfoil_type::control_point_type cst_airfoil_control_point_type;
    typedef eli::geom::curve::piecewise_cst_airfoil_fitter<data_type, 3, tolerance_type> airfoil_fitter_type;


    tolerance_type tol;

  private:
    void create_negative_x_airfoil_points(std::vector<point_type, Eigen::aligned_allocator<point_type>> &upt,
                                          std::vector<point_type, Eigen::aligned_allocator<point_type>> &lpt)
    {
      upt.resize(176);
      upt[ 0] << 0, 0, 0;
      upt[ 1] << -0.000422489, 0.00193697, 0; 
      upt[ 2] << -0.000445128, 0.00279966, 0; 
      upt[ 3] << -0.000407684, 0.00349485, 0; 
      upt[ 4] << -0.000336627, 0.00409544, 0; 
      upt[ 5] << -0.000241683, 0.00462858, 0; 
      upt[ 6] << -0.000133687, 0.00512366, 0; 
      upt[ 7] << -1.24082e-05, 0.00558071, 0; 
      upt[ 8] << 0.000117815, 0.00601158, 0; 
      upt[ 9] << 0.000255022, 0.00642174, 0; 
      upt[10] << 0.00039983, 0.00680984, 0; 
      upt[11] << 0.000549374, 0.00718361, 0; 
      upt[12] << 0.000703644, 0.00754324, 0; 
      upt[13] << 0.00086284, 0.00788835, 0; 
      upt[14] << 0.00102556, 0.00822275, 0; 
      upt[15] << 0.00119112, 0.00854833, 0; 
      upt[16] << 0.00136031, 0.00886311, 0; 
      upt[17] << 0.00153248, 0.00916889, 0; 
      upt[18] << 0.00170688, 0.00946767, 0; 
      upt[19] << 0.00188341, 0.00975979, 0; 
      upt[20] << 0.0020625, 0.0100442, 0; 
      upt[21] << 0.00227073, 0.0103321, 0; 
      upt[22] << 0.00248078, 0.0106142, 0; 
      upt[23] << 0.00269242, 0.0108914, 0; 
      upt[24] << 0.00290607, 0.0111622, 0; 
      upt[25] << 0.00312148, 0.0114275, 0; 
      upt[26] << 0.00333831, 0.011688, 0; 
      upt[27] << 0.00355623, 0.0119448, 0; 
      upt[28] << 0.00377518, 0.012198, 0; 
      upt[29] << 0.00399545, 0.0124468, 0; 
      upt[30] << 0.00421689, 0.0126915, 0; 
      upt[31] << 0.00443933, 0.0129327, 0; 
      upt[32] << 0.00466262, 0.0131709, 0; 
      upt[33] << 0.00488663, 0.0134062, 0; 
      upt[34] << 0.00511159, 0.0136383, 0; 
      upt[35] << 0.00533744, 0.013867, 0; 
      upt[36] << 0.00556409, 0.0140929, 0; 
      upt[37] << 0.00579144, 0.014316, 0; 
      upt[38] << 0.00601939, 0.0145368, 0; 
      upt[39] << 0.00624793, 0.0147553, 0; 
      upt[40] << 0.00647719, 0.0149711, 0; 
      upt[41] << 0.00670712, 0.0151843, 0; 
      upt[42] << 0.00693765, 0.0153951, 0; 
      upt[43] << 0.00716873, 0.0156038, 0; 
      upt[44] << 0.00740028, 0.0158105, 0; 
      upt[45] << 0.00763227, 0.0160153, 0; 
      upt[46] << 0.00786482, 0.0162179, 0; 
      upt[47] << 0.00809791, 0.0164183, 0; 
      upt[48] << 0.00833148, 0.0166167, 0; 
      upt[49] << 0.00856549, 0.0168133, 0; 
      upt[50] << 0.00879988, 0.0170082, 0; 
      upt[51] << 0.0111627, 0.0188729, 0; 
      upt[52] << 0.013554, 0.020605, 0; 
      upt[53] << 0.0159659, 0.0222327, 0; 
      upt[54] << 0.0183938, 0.0237752, 0; 
      upt[55] << 0.0208335, 0.0252484, 0; 
      upt[56] << 0.0232833, 0.0266608, 0; 
      upt[57] << 0.0257416, 0.0280199, 0; 
      upt[58] << 0.0282068, 0.0293327, 0; 
      upt[59] << 0.0306787, 0.030601, 0; 
      upt[60] << 0.0331563, 0.0318294, 0; 
      upt[61] << 0.0356389, 0.0330216, 0; 
      upt[62] << 0.0381262, 0.0341795, 0; 
      upt[63] << 0.0406179, 0.0353037, 0; 
      upt[64] << 0.0431139, 0.0363959, 0; 
      upt[65] << 0.0456134, 0.0374593, 0; 
      upt[66] << 0.0481163, 0.038495, 0; 
      upt[67] << 0.0506225, 0.039503, 0; 
      upt[68] << 0.0531317, 0.0404853, 0; 
      upt[69] << 0.0556437, 0.0414432, 0; 
      upt[70] << 0.0581581, 0.0423777, 0; 
      upt[71] << 0.0606748, 0.0432901, 0; 
      upt[72] << 0.0631937, 0.0441814, 0; 
      upt[73] << 0.0657145, 0.0450524, 0; 
      upt[74] << 0.0682371, 0.0459041, 0; 
      upt[75] << 0.0707614, 0.0467373, 0; 
      upt[76] << 0.0732873, 0.0475525, 0; 
      upt[77] << 0.0758147, 0.0483505, 0; 
      upt[78] << 0.0783433, 0.0491321, 0; 
      upt[79] << 0.0808731, 0.0498991, 0; 
      upt[80] << 0.0834041, 0.0506503, 0; 
      upt[81] << 0.0859363, 0.0513853, 0; 
      upt[82] << 0.0884696, 0.052106, 0; 
      upt[83] << 0.0910039, 0.0528126, 0; 
      upt[84] << 0.0935388, 0.053507, 0; 
      upt[85] << 0.0960746, 0.0541883, 0; 
      upt[86] << 0.106226, 0.0567834, 0; 
      upt[87] << 0.116388, 0.0591904, 0; 
      upt[88] << 0.126557, 0.0614272, 0; 
      upt[89] << 0.136732, 0.0635049, 0; 
      upt[90] << 0.146913, 0.0654357, 0; 
      upt[91] << 0.157096, 0.0672239, 0; 
      upt[92] << 0.167283, 0.0688861, 0; 
      upt[93] << 0.17747, 0.0704199, 0; 
      upt[94] << 0.187659, 0.0718321, 0; 
      upt[95] << 0.197848, 0.0731335, 0; 
      upt[96] << 0.208037, 0.0743252, 0; 
      upt[97] << 0.218224, 0.0754118, 0; 
      upt[98] << 0.22841, 0.0763955, 0; 
      upt[99] << 0.238594, 0.0772856, 0; 
      upt[100] << 0.248775, 0.0780791, 0; 
      upt[101] << 0.258954, 0.0787866, 0; 
      upt[102] << 0.269129, 0.0794063, 0; 
      upt[103] << 0.279301, 0.0799412, 0; 
      upt[104] << 0.28947, 0.0803908, 0; 
      upt[105] << 0.299634, 0.0807533, 0; 
      upt[106] << 0.309794, 0.0810379, 0; 
      upt[107] << 0.31995, 0.0812405, 0; 
      upt[108] << 0.330101, 0.0813606, 0; 
      upt[109] << 0.340247, 0.0814029, 0; 
      upt[110] << 0.350388, 0.081363, 0; 
      upt[111] << 0.360524, 0.08125, 0; 
      upt[112] << 0.370654, 0.0810606, 0; 
      upt[113] << 0.380778, 0.0807929, 0; 
      upt[114] << 0.390897, 0.0804575, 0; 
      upt[115] << 0.401009, 0.0800512, 0; 
      upt[116] << 0.411115, 0.0795743, 0; 
      upt[117] << 0.421215, 0.0790297, 0; 
      upt[118] << 0.431309, 0.0784217, 0; 
      upt[119] << 0.441397, 0.0777493, 0; 
      upt[120] << 0.451477, 0.0770136, 0; 
      upt[121] << 0.461552, 0.0762219, 0; 
      upt[122] << 0.47162, 0.0753726, 0; 
      upt[123] << 0.481682, 0.0744686, 0; 
      upt[124] << 0.491738, 0.0735123, 0; 
      upt[125] << 0.501787, 0.0725022, 0; 
      upt[126] << 0.51183, 0.0714447, 0; 
      upt[127] << 0.521867, 0.0703408, 0; 
      upt[128] << 0.531898, 0.0691878, 0; 
      upt[129] << 0.541923, 0.0679956, 0; 
      upt[130] << 0.551942, 0.0667602, 0; 
      upt[131] << 0.561956, 0.0654839, 0; 
      upt[132] << 0.571964, 0.0641705, 0; 
      upt[133] << 0.581966, 0.0628202, 0; 
      upt[134] << 0.591963, 0.0614357, 0; 
      upt[135] << 0.601955, 0.0600157, 0; 
      upt[136] << 0.611942, 0.0585651, 0; 
      upt[137] << 0.621924, 0.0570808, 0; 
      upt[138] << 0.631901, 0.0555703, 0; 
      upt[139] << 0.641874, 0.05403, 0; 
      upt[140] << 0.651842, 0.0524645, 0; 
      upt[141] << 0.661806, 0.0508726, 0; 
      upt[142] << 0.671766, 0.0492608, 0; 
      upt[143] << 0.681723, 0.0476291, 0; 
      upt[144] << 0.691676, 0.0459772, 0; 
      upt[145] << 0.701625, 0.0443116, 0; 
      upt[146] << 0.711572, 0.0426289, 0; 
      upt[147] << 0.721516, 0.0409378, 0; 
      upt[148] << 0.731457, 0.0392319, 0; 
      upt[149] << 0.741396, 0.0375199, 0; 
      upt[150] << 0.751333, 0.0358009, 0; 
      upt[151] << 0.761268, 0.0340798, 0; 
      upt[152] << 0.771201, 0.0323573, 0; 
      upt[153] << 0.781134, 0.0306363, 0; 
      upt[154] << 0.791066, 0.0289179, 0; 
      upt[155] << 0.800997, 0.0272052, 0; 
      upt[156] << 0.810927, 0.0255015, 0; 
      upt[157] << 0.820858, 0.0238117, 0; 
      upt[158] << 0.83079, 0.0221374, 0; 
      upt[159] << 0.840721, 0.0204814, 0; 
      upt[160] << 0.850654, 0.0188462, 0; 
      upt[161] << 0.860589, 0.0172362, 0; 
      upt[162] << 0.870525, 0.0156555, 0; 
      upt[163] << 0.880463, 0.0141026, 0; 
      upt[164] << 0.890403, 0.0125865, 0; 
      upt[165] << 0.900346, 0.0111084, 0; 
      upt[166] << 0.910291, 0.00967471, 0; 
      upt[167] << 0.92024, 0.0082886, 0; 
      upt[168] << 0.930193, 0.00695867, 0; 
      upt[169] << 0.940149, 0.00569217, 0; 
      upt[170] << 0.950111, 0.00449467, 0; 
      upt[171] << 0.960077, 0.00337381, 0; 
      upt[172] << 0.970048, 0.00234303, 0; 
      upt[173] << 0.980025, 0.00141442, 0; 
      upt[174] << 0.990008, 0.000612813, 0; 
      upt[175] << 0.999999, 1.70117e-07, 0;

      lpt.resize(176);
      lpt[ 0] << 0, 0, 0;
      lpt[ 1] << 0.000922489, -0.00158711, 0; 
      lpt[ 2] << 0.00144513, -0.0021552, 0; 
      lpt[ 3] << 0.00190768, -0.00257667, 0; 
      lpt[ 4] << 0.00233663, -0.00291713, 0; 
      lpt[ 5] << 0.00274168, -0.00320026, 0; 
      lpt[ 6] << 0.00313369, -0.00345342, 0; 
      lpt[ 7] << 0.00351241, -0.00367527, 0; 
      lpt[ 8] << 0.00388219, -0.00387671, 0; 
      lpt[ 9] << 0.00424498, -0.00406249, 0; 
      lpt[10] << 0.00460017, -0.00423069, 0; 
      lpt[11] << 0.00495063, -0.00438861, 0; 
      lpt[12] << 0.00529636, -0.00453608, 0; 
      lpt[13] << 0.00563716, -0.00467241, 0; 
      lpt[14] << 0.00597444, -0.00480113, 0; 
      lpt[15] << 0.00630888, -0.00492395, 0; 
      lpt[16] << 0.00663969, -0.00503867, 0; 
      lpt[17] << 0.00696752, -0.00514693, 0; 
      lpt[18] << 0.00729312, -0.00525059, 0; 
      lpt[19] << 0.00761659, -0.00534986, 0; 
      lpt[20] << 0.0079375, -0.00544354, 0; 
      lpt[21] << 0.00822927, -0.00554276, 0; 
      lpt[22] << 0.00851922, -0.00563824, 0; 
      lpt[23] << 0.00880759, -0.0057305, 0; 
      lpt[24] << 0.00909394, -0.00581828, 0; 
      lpt[25] << 0.00937852, -0.00590216, 0; 
      lpt[26] << 0.00966169, -0.005983, 0; 
      lpt[27] << 0.00994377, -0.00606164, 0; 
      lpt[28] << 0.0102248, -0.00613824, 0; 
      lpt[29] << 0.0105046, -0.00621186, 0; 
      lpt[30] << 0.0107831, -0.00628288, 0; 
      lpt[31] << 0.0110607, -0.00635173, 0; 
      lpt[32] << 0.0113374, -0.00641885, 0; 
      lpt[33] << 0.0116134, -0.00648451, 0; 
      lpt[34] << 0.0118884, -0.00654806, 0; 
      lpt[35] << 0.0121626, -0.00660958, 0; 
      lpt[36] << 0.0124359, -0.00666934, 0; 
      lpt[37] << 0.0127086, -0.0067276, 0; 
      lpt[38] << 0.0129806, -0.00678461, 0; 
      lpt[39] << 0.0132521, -0.00684041, 0; 
      lpt[40] << 0.0135228, -0.00689454, 0; 
      lpt[41] << 0.0137929, -0.00694715, 0; 
      lpt[42] << 0.0140623, -0.00699839, 0; 
      lpt[43] << 0.0143313, -0.00704846, 0; 
      lpt[44] << 0.0145997, -0.00709751, 0; 
      lpt[45] << 0.0148677, -0.00714567, 0; 
      lpt[46] << 0.0151352, -0.00719252, 0; 
      lpt[47] << 0.0154021, -0.00723811, 0; 
      lpt[48] << 0.0156685, -0.00728258, 0; 
      lpt[49] << 0.0159345, -0.00732609, 0; 
      lpt[50] << 0.0162001, -0.00736878, 0; 
      lpt[51] << 0.0188373, -0.00775519, 0; 
      lpt[52] << 0.021446, -0.00807966, 0; 
      lpt[53] << 0.024034, -0.00836086, 0; 
      lpt[54] << 0.0266062, -0.00861053, 0; 
      lpt[55] << 0.0291665, -0.00883899, 0; 
      lpt[56] << 0.0317167, -0.00905016, 0; 
      lpt[57] << 0.0342584, -0.00924785, 0; 
      lpt[58] << 0.0367932, -0.00943584, 0; 
      lpt[59] << 0.0393213, -0.00961339, 0; 
      lpt[60] << 0.0418437, -0.00978267, 0; 
      lpt[61] << 0.0443611, -0.00994561, 0; 
      lpt[62] << 0.0468738, -0.0101022, 0; 
      lpt[63] << 0.0493821, -0.0102518, 0; 
      lpt[64] << 0.0518861, -0.0103944, 0; 
      lpt[65] << 0.0543866, -0.0105322, 0; 
      lpt[66] << 0.0568837, -0.0106651, 0; 
      lpt[67] << 0.0593775, -0.0107921, 0; 
      lpt[68] << 0.0618683, -0.0109144, 0; 
      lpt[69] << 0.0643564, -0.0110323, 0; 
      lpt[70] << 0.0668419, -0.0111462, 0; 
      lpt[71] << 0.0693252, -0.0112565, 0; 
      lpt[72] << 0.0718063, -0.0113638, 0; 
      lpt[73] << 0.0742855, -0.0114681, 0; 
      lpt[74] << 0.0767629, -0.0115699, 0; 
      lpt[75] << 0.0792386, -0.0116694, 0; 
      lpt[76] << 0.0817127, -0.0117667, 0; 
      lpt[77] << 0.0841853, -0.0118621, 0; 
      lpt[78] << 0.0866567, -0.0119561, 0; 
      lpt[79] << 0.0891269, -0.0120498, 0; 
      lpt[80] << 0.0915959, -0.0121417, 0; 
      lpt[81] << 0.0940636, -0.0122313, 0; 
      lpt[82] << 0.0965303, -0.0123197, 0; 
      lpt[83] << 0.0989961, -0.0124073, 0; 
      lpt[84] << 0.101461, -0.0124953, 0; 
      lpt[85] << 0.103925, -0.0125826, 0; 
      lpt[86] << 0.113774, -0.0129204, 0; 
      lpt[87] << 0.123612, -0.0132483, 0; 
      lpt[88] << 0.133443, -0.0135705, 0; 
      lpt[89] << 0.143268, -0.0138864, 0; 
      lpt[90] << 0.153087, -0.0141981, 0; 
      lpt[91] << 0.162903, -0.014501, 0; 
      lpt[92] << 0.172717, -0.0148038, 0; 
      lpt[93] << 0.18253, -0.0150971, 0; 
      lpt[94] << 0.192341, -0.0153815, 0; 
      lpt[95] << 0.202152, -0.0156619, 0; 
      lpt[96] << 0.211963, -0.0159343, 0; 
      lpt[97] << 0.221776, -0.0161985, 0; 
      lpt[98] << 0.23159, -0.0164523, 0; 
      lpt[99] << 0.241406, -0.016701, 0; 
      lpt[100] << 0.251225, -0.0169378, 0; 
      lpt[101] << 0.261046, -0.0171697, 0; 
      lpt[102] << 0.270871, -0.0173916, 0; 
      lpt[103] << 0.280699, -0.0176036, 0; 
      lpt[104] << 0.29053, -0.0178023, 0; 
      lpt[105] << 0.300366, -0.0179828, 0; 
      lpt[106] << 0.310206, -0.0181521, 0; 
      lpt[107] << 0.32005, -0.0183035, 0; 
      lpt[108] << 0.329899, -0.0184341, 0; 
      lpt[109] << 0.339752, -0.0185464, 0; 
      lpt[110] << 0.349612, -0.0186338, 0; 
      lpt[111] << 0.359476, -0.0187035, 0; 
      lpt[112] << 0.369346, -0.0187502, 0; 
      lpt[113] << 0.379222, -0.01877, 0; 
      lpt[114] << 0.389103, -0.018772, 0; 
      lpt[115] << 0.398991, -0.0187508, 0; 
      lpt[116] << 0.408884, -0.0187053, 0; 
      lpt[117] << 0.418784, -0.0186367, 0; 
      lpt[118] << 0.428691, -0.0185479, 0; 
      lpt[119] << 0.438603, -0.0184361, 0; 
      lpt[120] << 0.448522, -0.018301, 0; 
      lpt[121] << 0.458448, -0.0181485, 0; 
      lpt[122] << 0.468379, -0.0179754, 0; 
      lpt[123] << 0.478317, -0.0177835, 0; 
      lpt[124] << 0.488262, -0.0175737, 0; 
      lpt[125] << 0.498212, -0.0173433, 0; 
      lpt[126] << 0.508169, -0.0170972, 0; 
      lpt[127] << 0.518132, -0.0168354, 0; 
      lpt[128] << 0.528101, -0.0165538, 0; 
      lpt[129] << 0.538076, -0.0162611, 0; 
      lpt[130] << 0.548057, -0.0159521, 0; 
      lpt[131] << 0.558044, -0.015628, 0; 
      lpt[132] << 0.568036, -0.0152914, 0; 
      lpt[133] << 0.578033, -0.0149413, 0; 
      lpt[134] << 0.588036, -0.0145793, 0; 
      lpt[135] << 0.598044, -0.0142029, 0; 
      lpt[136] << 0.608057, -0.0138161, 0; 
      lpt[137] << 0.618075, -0.0134145, 0; 
      lpt[138] << 0.628098, -0.0130046, 0; 
      lpt[139] << 0.638126, -0.0125817, 0; 
      lpt[140] << 0.648157, -0.0121494, 0; 
      lpt[141] << 0.658193, -0.0117052, 0; 
      lpt[142] << 0.668233, -0.0112546, 0; 
      lpt[143] << 0.678277, -0.0107968, 0; 
      lpt[144] << 0.688324, -0.01033, 0; 
      lpt[145] << 0.698374, -0.00985995, 0; 
      lpt[146] << 0.708428, -0.00938201, 0; 
      lpt[147] << 0.718484, -0.00890385, 0; 
      lpt[148] << 0.728543, -0.00841794, 0; 
      lpt[149] << 0.738604, -0.00793194, 0; 
      lpt[150] << 0.748667, -0.00744382, 0; 
      lpt[151] << 0.758731, -0.0069573, 0; 
      lpt[152] << 0.768798, -0.00647191, 0; 
      lpt[153] << 0.778865, -0.0059895, 0; 
      lpt[154] << 0.788934, -0.00550996, 0; 
      lpt[155] << 0.799002, -0.00503528, 0; 
      lpt[156] << 0.809072, -0.00456732, 0; 
      lpt[157] << 0.819141, -0.0041098, 0; 
      lpt[158] << 0.82921, -0.00366305, 0; 
      lpt[159] << 0.839278, -0.00322864, 0; 
      lpt[160] << 0.849345, -0.00280748, 0; 
      lpt[161] << 0.85941, -0.00240273, 0; 
      lpt[162] << 0.869474, -0.00201679, 0; 
      lpt[163] << 0.879536, -0.00164669, 0; 
      lpt[164] << 0.889596, -0.00129986, 0; 
      lpt[165] << 0.899653, -0.000975561, 0; 
      lpt[166] << 0.909708, -0.000678325, 0; 
      lpt[167] << 0.919759, -0.000409356, 0; 
      lpt[168] << 0.929806, -0.000174991, 0; 
      lpt[169] << 0.939849, 2.00381e-05, 0; 
      lpt[170] << 0.949888, 0.000172962, 0; 
      lpt[171] << 0.959922, 0.000279405, 0; 
      lpt[172] << 0.969951, 0.000329852, 0; 
      lpt[173] << 0.979974, 0.000317196, 0; 
      lpt[174] << 0.989991, 0.000223591, 0; 
      lpt[175] << 0.999999, -1.17942e-07, 0;
    }

    void create_airfoil_points(std::vector<point_type, Eigen::aligned_allocator<point_type>> &pts)
    {
      pts.resize(40);
      // lower surface
      pts[ 0] << 1.0000,-0.0020, 0;
      pts[ 1] << 0.9012,-0.0087, 0;
      pts[ 2] << 0.8078,-0.0247, 0;
      pts[ 3] << 0.7017,-0.0453, 0;
      pts[ 4] << 0.6178,-0.0572, 0;
      pts[ 5] << 0.4123,-0.0684, 0;
      pts[ 6] << 0.3545,-0.0678, 0;
      pts[ 7] << 0.2986,-0.0657, 0;
      pts[ 8] << 0.2453,-0.0621, 0;
      pts[ 9] << 0.1499,-0.0511, 0;
      pts[10] << 0.1029,-0.0441, 0;
      pts[11] << 0.0741,-0.0365, 0;
      pts[12] << 0.0451,-0.0284, 0;
      pts[13] << 0.0143,-0.0112, 0;
      pts[14] << 0.0076,-0.0116, 0;
      pts[15] << 0.0029,-0.0077, 0;
      pts[16] << 0.0000, 0.0000, 0;
      // upper surface
      pts[17] << 0.0013, 0.0030, 0;
      pts[18] << 0.0045, 0.0094, 0;
      pts[19] << 0.0096, 0.0138, 0;
      pts[20] << 0.0165, 0.0183, 0;
      pts[21] << 0.0252, 0.0228, 0;
      pts[22] << 0.0477, 0.0315, 0;
      pts[23] << 0.0767, 0.0397, 0;
      pts[24] << 0.1118, 0.0473, 0;
      pts[25] << 0.1524, 0.0539, 0;
      pts[26] << 0.1980, 0.0593, 0;
      pts[27] << 0.2979, 0.0636, 0;
      pts[28] << 0.3015, 0.0665, 0;
      pts[29] << 0.3578, 0.0680, 0;
      pts[30] << 0.4160, 0.0680, 0;
      pts[31] << 0.4455, 0.0675, 0;
      pts[32] << 0.5049, 0.0652, 0;
      pts[33] << 0.5930, 0.0585, 0;
      pts[34] << 0.6501, 0.0514, 0;
      pts[35] << 0.7050, 0.0416, 0;
      pts[36] << 0.7623, 0.0297, 0;
      pts[37] << 0.8168, 0.0221, 0;
      pts[38] << 0.9074, 0.0108, 0;
      pts[39] << 1.0000, 0.0050, 0;
    }

    void create_airfoil_points(std::vector<point_type, Eigen::aligned_allocator<point_type>> &upt,
                               std::vector<point_type, Eigen::aligned_allocator<point_type>> &lpt)
    {
      // upper surface
      upt.resize(24);
      upt[ 0] << 0.0000, 0.0000, 0;
      upt[ 1] << 0.0013, 0.0030, 0;
      upt[ 2] << 0.0045, 0.0094, 0;
      upt[ 3] << 0.0096, 0.0138, 0;
      upt[ 4] << 0.0165, 0.0183, 0;
      upt[ 5] << 0.0252, 0.0228, 0;
      upt[ 6] << 0.0477, 0.0315, 0;
      upt[ 7] << 0.0767, 0.0397, 0;
      upt[ 8] << 0.1118, 0.0473, 0;
      upt[ 9] << 0.1524, 0.0539, 0;
      upt[10] << 0.1980, 0.0593, 0;
      upt[11] << 0.2979, 0.0636, 0;
      upt[12] << 0.3015, 0.0665, 0;
      upt[13] << 0.3578, 0.0680, 0;
      upt[14] << 0.4160, 0.0680, 0;
      upt[15] << 0.4455, 0.0675, 0;
      upt[16] << 0.5049, 0.0652, 0;
      upt[17] << 0.5930, 0.0585, 0;
      upt[18] << 0.6501, 0.0514, 0;
      upt[19] << 0.7050, 0.0416, 0;
      upt[20] << 0.7623, 0.0297, 0;
      upt[21] << 0.8168, 0.0221, 0;
      upt[22] << 0.9074, 0.0108, 0;
      upt[23] << 1.0000, 0.0050, 0;

      // lower surface
      lpt.resize(17);
      lpt[ 0] << 0.0000, 0.0000, 0;
      lpt[ 1] << 0.0029,-0.0077, 0;
      lpt[ 2] << 0.0076,-0.0116, 0;
      lpt[ 3] << 0.0143,-0.0112, 0;
      lpt[ 4] << 0.0451,-0.0284, 0;
      lpt[ 5] << 0.0741,-0.0365, 0;
      lpt[ 6] << 0.1029,-0.0441, 0;
      lpt[ 7] << 0.1499,-0.0511, 0;
      lpt[ 8] << 0.2453,-0.0621, 0;
      lpt[ 9] << 0.2986,-0.0657, 0;
      lpt[10] << 0.3545,-0.0678, 0;
      lpt[11] << 0.4123,-0.0684, 0;
      lpt[12] << 0.6178,-0.0572, 0;
      lpt[13] << 0.7017,-0.0453, 0;
      lpt[14] << 0.8078,-0.0247, 0;
      lpt[15] << 0.9012,-0.0087, 0;
      lpt[16] << 1.0000,-0.0020, 0;
    }


    void create_airfoil_points(std::vector<point_type, Eigen::aligned_allocator<point_type>> &upt,
                               std::vector<point_type, Eigen::aligned_allocator<point_type>> &lpt,
                               point_type & le_pt, point_type &te_pt)
    {
      // get the airfoil points
      create_airfoil_points(upt, lpt);

      // adjust the lower surface points to put trailing edge half-way
      // between upper and lower trailing edge points
      data_type y_offset(-0.003);

      for (index_type i=0; i<static_cast<index_type>(lpt.size()); ++i)
      {
        lpt[i].y()+=lpt[i].x()*y_offset;
      }

      // transform airfoil points
      data_type scale(2.5);
      data_type theta(30*eli::constants::math<data_type>::pi()/180);
      le_pt << 1, 2, 0;
      te_pt << le_pt.x()+std::cos(theta), le_pt.y()+std::sin(theta), 0;

      for (index_type i=0; i<static_cast<index_type>(upt.size()); ++i)
      {
        // (1) scale
        upt[i] = scale*upt[i];

        // (2) rotate
        data_type ptx(upt[i].x()), pty(upt[i].y()), ptz(upt[i].z()), ct(std::cos(theta)), st(std::sin(theta));
        upt[i] << (ptx*ct-pty*st), (ptx*st+pty*ct), ptz;

        // (3) translate
        upt[i] = upt[i]+le_pt;
      }
      for (index_type i=0; i<static_cast<index_type>(lpt.size()); ++i)
      {
        // (1) scale
        lpt[i] = scale*lpt[i];

        // (2) rotate
        data_type ptx(lpt[i].x()), pty(lpt[i].y()), ptz(lpt[i].z()), ct(std::cos(theta)), st(std::sin(theta));
        lpt[i] << (ptx*ct-pty*st), (ptx*st+pty*ct), ptz;

        // (3) translate
        lpt[i] = lpt[i]+le_pt;
      }
    }

    void create_cst_airfoil_points(std::vector<point_type, Eigen::aligned_allocator<point_type>> &upt,
                                   std::vector<point_type, Eigen::aligned_allocator<point_type>> &lpt)
    {
      std::vector<cst_airfoil_control_point_type, Eigen::aligned_allocator<cst_airfoil_control_point_type>> cpu, cpl;

      create_cst_airfoil_points(upt, lpt, cpu, cpl);
    }

    void create_cst_airfoil_points(std::vector<point_type, Eigen::aligned_allocator<point_type>> &upt,
                                   std::vector<point_type, Eigen::aligned_allocator<point_type>> &lpt,
                                   std::vector<cst_airfoil_control_point_type, Eigen::aligned_allocator<cst_airfoil_control_point_type>> &cpu,
                                   std::vector<cst_airfoil_control_point_type, Eigen::aligned_allocator<cst_airfoil_control_point_type>> &cpl)
    {
      // create an asymmetric CST airfoil
      typedef typename cst_airfoil_type::point_type cst_airfoil_point_type;
      typedef typename cst_airfoil_type::data_type cst_airfoil_data_type;
      typedef typename cst_airfoil_type::index_type cst_airfoil_index_type;

      cst_airfoil_type cst(7, 5);
      cst_airfoil_data_type dteu(0.007), dtel(0.005);
      cst_airfoil_index_type i;

      // set the control points
      cpu.resize(8);
      cpu[0] << static_cast<data_type>( 0.17);
      cpu[1] << static_cast<data_type>( 0.16);
      cpu[2] << static_cast<data_type>( 0.16);
      cpu[3] << static_cast<data_type>( 0.14);
      cpu[4] << static_cast<data_type>( 0.15);
      cpu[5] << static_cast<data_type>( 0.14);
      cpu[6] << static_cast<data_type>( 0.14);
      cpu[7] << static_cast<data_type>( 0.14);
      cpl.resize(6);
      cpl[0] << static_cast<data_type>(-0.17);
      cpl[1] << static_cast<data_type>( 0.08);
      cpl[2] << static_cast<data_type>( 0.06);
      cpl[3] << static_cast<data_type>( 0.00);
      cpl[4] << static_cast<data_type>( 0.00);
      cpl[5] << static_cast<data_type>( 0.05);
      for (i=0; i<=cst.upper_degree(); ++i)
      {
        cst.set_upper_control_point(cpu[i], i);
      }
      for (i=0; i<=cst.lower_degree(); ++i)
      {
        cst.set_lower_control_point(cpl[i], i);
      }

      // set the trailing edge thickness of CST airfoil
      cst.set_trailing_edge_thickness(dteu, dtel);

      // sample the lower and upper surfaces
      upt.resize(24);
      for (i=0; i<static_cast<cst_airfoil_index_type>(upt.size()); ++i)
      {
        cst_airfoil_point_type pt(cst.f(static_cast<data_type>(i)/(upt.size()-1)));
        upt[i] << pt.x(), pt.y(), 0;
//        std::cout << "upt[" << i << "]=" << upt[i] << std::endl;
      }
      lpt.resize(20);
      for (i=0; i<static_cast<cst_airfoil_index_type>(lpt.size()); ++i)
      {
        cst_airfoil_point_type pt(cst.f(-static_cast<data_type>(i)/(lpt.size()-1)));
        lpt[i] << pt.x(), pt.y(), 0;
//        std::cout << "lpt[" << i << "]=" << lpt[i] << std::endl;
      }
    }

  protected:
    void AddTests(const float &)
    {
      // add the tests
//      std::cout << "%% remove comments here" << std::endl;
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<float>::create_airfoil_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<float>::fit_airfoil_to_cst_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<float>::fit_airfoil_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<float>::fit_negative_x_airfoil_test);
    }
    void AddTests(const double &)
    {
      // add the tests
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<double>::create_airfoil_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<double>::fit_airfoil_to_cst_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<double>::fit_airfoil_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<double>::fit_negative_x_airfoil_test);
    }
    void AddTests(const long double &)
    {
      // add the tests
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<long double>::create_airfoil_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<long double>::fit_airfoil_to_cst_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<long double>::fit_airfoil_test);
      TEST_ADD(piecewise_cst_airfoil_creator_test_suite<long double>::fit_negative_x_airfoil_test);
    }

  public:
    piecewise_cst_airfoil_creator_test_suite() : tol()
    {
      AddTests(data__());
    }
    ~piecewise_cst_airfoil_creator_test_suite()
    {
    }

  private:

    void create_airfoil_test()
    {
      typedef eli::geom::curve::piecewise_cst_airfoil_creator<data__, 3, tolerance_type> airfoil_creator_type;

      airfoil_creator_type pcst;
      piecewise_curve_type pc;
      cst_airfoil_type cst(7);
      cst_airfoil_control_point_type cp[8];
      data_type dte(2*0.00126), ti, t0, t1, t2, t[6];
      point_type pt_out, pt_ref;
      typename cst_airfoil_type::point_type pt2_ref;
      index_type i;
      bool rtn_flag;

      // set the control points
      cp[0] << static_cast<data_type>(0.170987592880629);
      cp[1] << static_cast<data_type>(0.157286894410384);
      cp[2] << static_cast<data_type>(0.162311658384540);
      cp[3] << static_cast<data_type>(0.143623187913493);
      cp[4] << static_cast<data_type>(0.149218456400780);
      cp[5] << static_cast<data_type>(0.137218405082418);
      cp[6] << static_cast<data_type>(0.140720628655908);
      cp[7] << static_cast<data_type>(0.141104769355436);
      for (i=0; i<=cst.upper_degree(); ++i)
      {
        cst.set_upper_control_point(cp[i], i);
        cst.set_lower_control_point(-cp[i], i);
      }

      // set the trailing edge thickness of CST airfoil
      cst.set_trailing_edge_thickness(dte);

      // set the parameterization
      t0=-1;
      t1=0;
      t2=1;

      // set the parameters to evaluate the tests
      t[0] = t0+(t2-t0)*static_cast<data_type>(0);
      t[1] = t0+(t2-t0)*static_cast<data_type>(0.1);
      t[2] = t0+(t2-t0)*static_cast<data_type>(0.27);
      t[3] = t0+(t2-t0)*static_cast<data_type>(0.5);
      t[4] = t0+(t2-t0)*static_cast<data_type>(0.73);
      t[5] = t0+(t2-t0)*static_cast<data_type>(1);

      // create curve
      rtn_flag=pcst.set_conditions(cst);
      TEST_ASSERT(rtn_flag);
      pcst.set_t0(t0);
      pcst.set_segment_dt(t1-t0, 0);
      pcst.set_segment_dt(t2-t1, 1);
      rtn_flag=pcst.create(pc);
      TEST_ASSERT(rtn_flag);

      // evaluate the points (note need to transform parameterization to match points
      i=0;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=1;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=2;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=3;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=4;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));
      i=5;
      ti=t[i];
      pt_out=pc.f((ti<0)?-std::sqrt(-ti) : std::sqrt(ti));
      ti=2*(t[i]-t0)/(t2-t0)-1;
      pt2_ref=cst.f(ti);
      pt_ref << pt2_ref(0), pt2_ref(1), 0;
      TEST_ASSERT(tol.approximately_equal(pt_out, pt_ref));

//      if (typeid(data_type)==typeid(float))
//      {
//        std::cout.flush();
//        eli::test::octave_start(1);
//        eli::test::octave_print(1, cst, "cst");
//        eli::test::octave_print(1, pc, "piecewise");
//        eli::test::octave_finish(1, false);
//      }
    }

    void fit_airfoil_to_cst_test()
    {
      {
        airfoil_fitter_type pcaf;
        cst_airfoil_type cst;
        std::vector<point_type, Eigen::aligned_allocator<point_type>> upt, lpt;
        std::vector<cst_airfoil_control_point_type, Eigen::aligned_allocator<cst_airfoil_control_point_type>> cpu_ref, cpl_ref;
        point_type pt, pt_ref;
        data_type t0, t1, t2;
        index_type degu, degl;
        bool rtn_flag;

        // get airfoil points
        create_cst_airfoil_points(upt, lpt, cpu_ref, cpl_ref);

        // set the parameterization
        t0=-1;
        t1=0;
        t2=1;

        // create curve
        degu=7;
        degl=5;
        rtn_flag=pcaf.set_conditions(upt.begin(), static_cast<index_type>(upt.size()), degu,
                                     lpt.begin(), static_cast<index_type>(lpt.size()), degl, false);
        TEST_ASSERT(rtn_flag);
        pcaf.set_t0(t0);
        pcaf.set_segment_dt(t1-t0, 0);
        pcaf.set_segment_dt(t2-t1, 1);


        Eigen::Matrix<data_type, 3, 3> transform_out;
        point_type translate_out;
        data_type actual_leading_edge_t;

        rtn_flag=pcaf.create(cst, transform_out, translate_out, actual_leading_edge_t);
        TEST_ASSERT(rtn_flag);

        // make sure got back correct angle and scale factor
        pt_ref.setZero();
        TEST_ASSERT(translate_out==pt_ref);
        TEST_ASSERT(tol.approximately_equal(transform_out, Eigen::Matrix<data_type, 3, 3>::Identity()));
        TEST_ASSERT(tol.approximately_equal(actual_leading_edge_t, 0));

        // cycle through CST control points to compare to the reference values
        cst_airfoil_control_point_type cp;
        index_type i;

        for (i=0; i<=cst.upper_degree(); ++i)
        {
          cp = cst.get_upper_control_point(i);
          std::string str("Error for Upper Control Point "+std::to_string(i));
          TEST_ASSERT_MSG(tol.approximately_equal(cp, cpu_ref[i]), str.data());
//          std::cout << "diff=" << (cp-cpu_ref[i]).norm() << std::endl;
        }

        for (i=0; i<=cst.lower_degree(); ++i)
        {
          cp = cst.get_lower_control_point(i);
          std::string str("Error for Lower Control Point "+std::to_string(i));
          TEST_ASSERT_MSG(tol.approximately_equal(cp, cpl_ref[i]), str.data());
//          std::cout << "diff=" << (cp-cpu_ref[i]).norm() << std::endl;
        }

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          // print out the upper surface points
//          for (size_t n=0; n<upt.size(); ++n)
//          {
//            std::string name("upt"); name+=std::to_string(n);
//            eli::test::octave_print(1, upt[n], name);
//          }
//          // print out the lower surface points
//          for (size_t n=0; n<lpt.size(); ++n)
//          {
//            std::string name("lpt"); name+=std::to_string(n);
//            eli::test::octave_print(1, lpt[n], name);
//          }
//          eli::test::octave_print(1, cst, "cst");
//          eli::test::octave_finish(1, true);
//        }
      }

      // fit to a known CST airfoil shape
      {
        airfoil_fitter_type pcaf;
        piecewise_curve_type pc;
        std::vector<point_type, Eigen::aligned_allocator<point_type>> upt, lpt;
        point_type pt, pt_ref;
        data_type t0, t1, t2;
        index_type degu, degl;
        bool rtn_flag;

        // get airfoil points
        create_cst_airfoil_points(upt, lpt);

        // set the parameterization
        t0=-1;
        t1=0;
        t2=1;

        // create curve
        degu=7;
        degl=5;
        rtn_flag=pcaf.set_conditions(upt.begin(), static_cast<index_type>(upt.size()), degu,
                                     lpt.begin(), static_cast<index_type>(lpt.size()), degl, false);
        TEST_ASSERT(rtn_flag);
        pcaf.set_t0(t0);
        pcaf.set_segment_dt(t1-t0, 0);
        pcaf.set_segment_dt(t2-t1, 1);
        rtn_flag=pcaf.create(pc);
        TEST_ASSERT(rtn_flag);

        // cycle through segments to get the control points to compare to the reference values
        control_point_type cp[32], cp_ref[32];
        curve_type crv;
        index_type i, ii;

        cp_ref[ 0] << static_cast<data_type>(1.00000000000000000), static_cast<data_type>(-0.00500000000000000), 0;
        cp_ref[ 1] << static_cast<data_type>(0.84615384615384600), static_cast<data_type>( 0.00346153846153849), 0;
        cp_ref[ 2] << static_cast<data_type>(0.70512820512820500), static_cast<data_type>(-0.00288461538461537), 0;
        cp_ref[ 3] << static_cast<data_type>(0.57692307692307700), static_cast<data_type>(-0.00288461538461520), 0;
        cp_ref[ 4] << static_cast<data_type>(0.46153846153846200), static_cast<data_type>( 0.01111888111888110), 0;
        cp_ref[ 5] << static_cast<data_type>(0.35897435897435900), static_cast<data_type>( 0.02306915306915300), 0;
        cp_ref[ 6] << static_cast<data_type>(0.26923076923076900), static_cast<data_type>( 0.01935314685314690), 0;
        cp_ref[ 7] << static_cast<data_type>(0.19230769230769200), static_cast<data_type>( 0.00146270396270403), 0;
        cp_ref[ 8] << static_cast<data_type>(0.12820512820512800), static_cast<data_type>(-0.01944444444444440), 0;
        cp_ref[ 9] << static_cast<data_type>(0.07692307692307690), static_cast<data_type>(-0.03283216783216790), 0;
        cp_ref[10] << static_cast<data_type>(0.03846153846153850), static_cast<data_type>(-0.03445804195804200), 0;
        cp_ref[11] << static_cast<data_type>(0.01282051282051280), static_cast<data_type>(-0.02621794871794880), 0;
        cp_ref[12] << static_cast<data_type>(0.00000000000000000), static_cast<data_type>(-0.01307692307692310), 0;
        cp_ref[13] << static_cast<data_type>(0.00000000000000000), static_cast<data_type>( 0.00000000000000000), 0;
        cp_ref[14] << static_cast<data_type>(0.00000000000000000), static_cast<data_type>( 0.00000000000000000), 0;
        cp_ref[15] << static_cast<data_type>(0.00000000000000000), static_cast<data_type>( 0.01000000000000000), 0;
        cp_ref[16] << static_cast<data_type>(0.00735294117647059), static_cast<data_type>( 0.02005147058823530), 0;
        cp_ref[17] << static_cast<data_type>(0.02205882352941180), static_cast<data_type>( 0.02980147058823530), 0;
        cp_ref[18] << static_cast<data_type>(0.04411764705882350), static_cast<data_type>( 0.03889705882352940), 0;
        cp_ref[19] << static_cast<data_type>(0.07352941176470590), static_cast<data_type>( 0.04703054298642530), 0;
        cp_ref[20] << static_cast<data_type>(0.11029411764705900), static_cast<data_type>( 0.05398472850678730), 0;
        cp_ref[21] << static_cast<data_type>(0.15441176470588200), static_cast<data_type>( 0.05961337926779100), 0;
        cp_ref[22] << static_cast<data_type>(0.20588235294117600), static_cast<data_type>( 0.06369210201563150), 0;
        cp_ref[23] << static_cast<data_type>(0.26470588235294100), static_cast<data_type>( 0.06573323735088460), 0;
        cp_ref[24] << static_cast<data_type>(0.33088235294117600), static_cast<data_type>( 0.06517508227067060), 0;
        cp_ref[25] << static_cast<data_type>(0.40441176470588200), static_cast<data_type>( 0.06229920814479620), 0;
        cp_ref[26] << static_cast<data_type>(0.48529411764705900), static_cast<data_type>( 0.05873642533936640), 0;
        cp_ref[27] << static_cast<data_type>(0.57352941176470600), static_cast<data_type>( 0.05454411764705930), 0;
        cp_ref[28] << static_cast<data_type>(0.66911764705882300), static_cast<data_type>( 0.04503676470588240), 0;
        cp_ref[29] << static_cast<data_type>(0.77205882352941200), static_cast<data_type>( 0.03525735294117730), 0;
        cp_ref[30] << static_cast<data_type>(0.88235294117647100), static_cast<data_type>( 0.02264705882353050), 0;
        cp_ref[31] << static_cast<data_type>(1.00000000000000000), static_cast<data_type>( 0.00700000000000000), 0;

        pc.get(crv, 0);
        ii=0;
        for (i=0; i<=crv.degree(); ++i, ++ii)
        {
          cp[ii]=crv.get_control_point(i);
        }
        pc.get(crv, 1);
        for (i=0; i<=crv.degree(); ++i, ++ii)
        {
          cp[ii]=crv.get_control_point(i);
        }

        for (i=0; i<32; ++i)
        {
          std::string str("Error for Control Point "+std::to_string(i));
          TEST_ASSERT_MSG(tol.approximately_equal(cp[i], cp_ref[i]), str.data());
//          std::cout << "cp_ref[" << std::setw(2) << i << "] << " << std::setprecision(15)
//                    << cp[i].x() << ", " << cp[i].y() << ", " << cp[i].z() << ";" << std::endl;
        }

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          // print out the upper surface points
//          for (size_t n=0; n<upt.size(); ++n)
//          {
//            std::string name("upt"); name+=std::to_string(n);
//            eli::test::octave_print(1, upt[n], name);
//          }
//          // print out the lower surface points
//          for (size_t n=0; n<lpt.size(); ++n)
//          {
//            std::string name("lpt"); name+=std::to_string(n);
//            eli::test::octave_print(1, lpt[n], name);
//          }
//          eli::test::octave_print(1, pc, "piecewise");
//          eli::test::octave_finish(1, true);
//        }
      }
    }

    void fit_airfoil_test()
    {
      // fit to simple two vector specification
      {
        airfoil_fitter_type pcaf;
        piecewise_curve_type pc;
        std::vector<point_type, Eigen::aligned_allocator<point_type>> upt, lpt;
        point_type pt, pt_ref;
        data_type t, t0, t1, t2;
        index_type degu, degl;
        bool rtn_flag;

        // get airfoil points
        create_airfoil_points(upt, lpt);

        // set the parameterization
        t0=-1;
        t1=0;
        t2=1;

        // create curve
        degu=8;
        degl=8;
        rtn_flag=pcaf.set_conditions(upt.begin(), static_cast<index_type>(upt.size()), degu,
                                     lpt.begin(), static_cast<index_type>(lpt.size()), degl, true);
        TEST_ASSERT(rtn_flag);
        pcaf.set_t0(t0);
        pcaf.set_segment_dt(t1-t0, 0);
        pcaf.set_segment_dt(t2-t1, 1);
        rtn_flag=pcaf.create(pc);
        TEST_ASSERT(rtn_flag);

        // test various points
        t  =-0.8;
        pt = pc.f(t);
        pt_ref << t*t, -0.05448192672, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);
        t  =-0.1;
        pt = pc.f(t);
        pt_ref << t*t, -0.0121366373, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);
        t  = 0.1;
        pt = pc.f(t);
        pt_ref << t*t, 0.01311318449, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);
        t  = 0.95;
        pt = pc.f(t);
        pt_ref << t*t, 0.01116413716, 0;
        TEST_ASSERT((pt-pt_ref).norm()<2e-4);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          // print out the upper surface points
//          for (size_t n=0; n<upt.size(); ++n)
//          {
//            std::string name("upt"); name+=std::to_string(n);
//            eli::test::octave_print(1, upt[n], name);
//          }
//          // print out the lower surface points
//          for (size_t n=0; n<lpt.size(); ++n)
//          {
//            std::string name("lpt"); name+=std::to_string(n);
//            eli::test::octave_print(1, lpt[n], name);
//          }
//          eli::test::octave_print(1, pc, "piecewise");
//          eli::test::octave_finish(1, true);
//        }
      }

      // scale, translate, and rotate points
      {
        airfoil_fitter_type pcaf;
        piecewise_curve_type pc;
        std::vector<point_type, Eigen::aligned_allocator<point_type>> upt, lpt;
        point_type pt, pt_ref, lept, tept;
        data_type t, t0, t1, t2;
        index_type degu, degl;
        bool rtn_flag;

        // get airfoil points
        create_airfoil_points(upt, lpt, lept, tept);

        // set the parameterization
        t0=-1;
        t1=0;
        t2=1;

        // create curve
        degu=8;
        degl=8;
        rtn_flag=pcaf.set_conditions(upt.begin(), static_cast<index_type>(upt.size()), degu,
                                     lpt.begin(), static_cast<index_type>(lpt.size()), degl, true);
        TEST_ASSERT(rtn_flag);
        pcaf.set_t0(t0);
        pcaf.set_segment_dt(t1-t0, 0);
        pcaf.set_segment_dt(t2-t1, 1);
        rtn_flag=pcaf.create(pc);
        TEST_ASSERT(rtn_flag);


        // test various points
        t  =-0.8;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 2.45614305445983, 2.67788624658704, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);

        t  =-0.1;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 1.03685893171540, 1.98615845755622, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);

        t  = 0.1;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 1.00525915448552, 2.04089087722623, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);

        t  = 0.95;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 2.93983922906606, 3.15259989674235, 0;
        TEST_ASSERT((pt-pt_ref).norm()<2e-4);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          // print out the upper surface points
//          for (size_t n=0; n<upt.size(); ++n)
//          {
//            std::string name("upt"); name+=std::to_string(n);
//            eli::test::octave_print(1, upt[n], name);
//          }
//          // print out the lower surface points
//          for (size_t n=0; n<lpt.size(); ++n)
//          {
//            std::string name("lpt"); name+=std::to_string(n);
//            eli::test::octave_print(1, lpt[n], name);
//          }
//          eli::test::octave_print(1, pc, "piecewise");
//          eli::test::octave_finish(1, true);
//        }
      }
    }

    void fit_negative_x_airfoil_test()
    {
      // fit to simple two vector specification
      {
        airfoil_fitter_type pcaf;
        piecewise_curve_type pc;
        std::vector<point_type, Eigen::aligned_allocator<point_type>> upt, lpt;
        point_type pt, pt_ref;
        data_type t, t0, t1, t2;
        index_type degu, degl;
        bool rtn_flag;

        // get airfoil points
        create_negative_x_airfoil_points(upt, lpt);

        // set the parameterization
        t0=-1;
        t1=0;
        t2=1;

        // create curve
        degu=8;
        degl=8;
        rtn_flag=pcaf.set_conditions(upt.begin(), static_cast<index_type>(upt.size()), degu,
                                     lpt.begin(), static_cast<index_type>(lpt.size()), degl, true);
        TEST_ASSERT(rtn_flag);
        pcaf.set_t0(t0);
        pcaf.set_segment_dt(t1-t0, 0);
        pcaf.set_segment_dt(t2-t1, 1);
        rtn_flag=pcaf.create(pc);
        TEST_ASSERT(rtn_flag);

        // test various points
        if (typeid(data_type)!=typeid(float)) // round off error is bad for float
        {
          t  =-0.8;
          pt = pc.f(t);
          pt_ref << 0.6466341658, -0.01154890063, 0;
//          std::cout << "pt_ref << " << std::setprecision(10) << pt.x() << " << " << std::setprecision(10) << pt.y() << ", 0;" << std::endl;
//          std::cout << "pt=" << pt << "\tpt_ref=" << pt_ref << "\tnorm=" << (pt-pt_ref).norm() << std::endl;
          TEST_ASSERT((pt-pt_ref).norm()<1e-5);
          t  =-0.1;
          pt = pc.f(t);
          pt_ref << 0.0137273553, -0.008180440922, 0;
//          std::cout << "pt_ref << " << std::setprecision(10) << pt.x() << " << " << std::setprecision(10) << pt.y() << ", 0;" << std::endl;
//          std::cout << "pt=" << pt << "\tpt_ref=" << pt_ref << "\tnorm=" << (pt-pt_ref).norm() << std::endl;
          TEST_ASSERT((pt-pt_ref).norm()<1e-5);
          t  = 0.1;
          pt = pc.f(t);
          pt_ref << 0.006058928226, 0.01409424987, 0;
//          std::cout << "pt_ref << " << std::setprecision(10) << pt.x() << " << " << std::setprecision(10) << pt.y() << ", 0;" << std::endl;
//          std::cout << "pt=" << pt << "\tpt_ref=" << pt_ref << "\tnorm=" << (pt-pt_ref).norm() << std::endl;
          TEST_ASSERT((pt-pt_ref).norm()<1e-5);
          t  = 0.95;
          pt = pc.f(t);
          pt_ref << 0.9004204782, 0.01095749739, 0;
//          std::cout << "pt_ref << " << std::setprecision(10) << pt.x() << " << " << std::setprecision(10) << pt.y() << ", 0;" << std::endl;
//          std::cout << "pt=" << pt << "\tpt_ref=" << pt_ref << "\tnorm=" << (pt-pt_ref).norm() << std::endl;
          TEST_ASSERT((pt-pt_ref).norm()<2e-4);

//        if (typeid(data_type)==typeid(double))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          // print out the upper surface points
//          for (size_t n=0; n<upt.size(); ++n)
//          {
//            std::string name("upt"); name+=std::to_string(n);
//            eli::test::octave_print(1, upt[n], name);
//          }
//          // print out the lower surface points
//          for (size_t n=0; n<lpt.size(); ++n)
//          {
//            std::string name("lpt"); name+=std::to_string(n);
//            eli::test::octave_print(1, lpt[n], name);
//          }
//          eli::test::octave_print(1, pc.f(0), "leading_edge");
//          eli::test::octave_print(1, pc, "piecewise");
//          eli::test::octave_finish(1, true);
//        }
        }
      }

#if 0
      // scale, translate, and rotate points
      {
        airfoil_fitter_type pcaf;
        piecewise_curve_type pc;
        std::vector<point_type, Eigen::aligned_allocator<point_type>> upt, lpt;
        point_type pt, pt_ref, lept, tept;
        data_type t, t0, t1, t2;
        index_type degu, degl;
        bool rtn_flag;

        // get airfoil points
        create_airfoil_points(upt, lpt, lept, tept);

        // set the parameterization
        t0=-1;
        t1=0;
        t2=1;

        // create curve
        degu=8;
        degl=8;
        rtn_flag=pcaf.set_conditions(upt.begin(), static_cast<index_type>(upt.size()), degu,
                                     lpt.begin(), static_cast<index_type>(lpt.size()), degl, true);
        TEST_ASSERT(rtn_flag);
        pcaf.set_t0(t0);
        pcaf.set_segment_dt(t1-t0, 0);
        pcaf.set_segment_dt(t2-t1, 1);
        rtn_flag=pcaf.create(pc);
        TEST_ASSERT(rtn_flag);


        // test various points
        t  =-0.8;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 2.45614305445983, 2.67788624658704, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);

        t  =-0.1;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 1.03685893171540, 1.98615845755622, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);

        t  = 0.1;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 1.00525915448552, 2.04089087722623, 0;
        TEST_ASSERT((pt-pt_ref).norm()<1e-5);

        t  = 0.95;
        pt = pc.f(t);
//        std::cout << "pt_ref << " << std::setprecision(15) << pt.x() << ", " << pt.y() << ", " << pt.z() << ";" << std::endl;
        pt_ref << 2.93983922906606, 3.15259989674235, 0;
        TEST_ASSERT((pt-pt_ref).norm()<2e-4);

//        if (typeid(data_type)==typeid(float))
//        {
//          std::cout.flush();
//          eli::test::octave_start(1);
//          // print out the upper surface points
//          for (size_t n=0; n<upt.size(); ++n)
//          {
//            std::string name("upt"); name+=std::to_string(n);
//            eli::test::octave_print(1, upt[n], name);
//          }
//          // print out the lower surface points
//          for (size_t n=0; n<lpt.size(); ++n)
//          {
//            std::string name("lpt"); name+=std::to_string(n);
//            eli::test::octave_print(1, lpt[n], name);
//          }
//          eli::test::octave_print(1, pc, "piecewise");
//          eli::test::octave_finish(1, true);
//        }
      }
#endif
    }
};

#endif
