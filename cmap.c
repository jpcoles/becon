#include <assert.h>
#include "cmap.h"

void cmap_astro(double v, int *r, int *g, int *b)
{
    static int ramp[256][3] =  {
  {0,0,255}, {0,0,254}, {1,0,253}, {2,0,252}, {3,0,250}, {4,0,249}, {5,0,248}, {6,0,247}, {7,0,247},
  {8,0,246}, {9,0,245}, {10,0,244}, {11,0,242}, {12,0,241}, {13,0,240}, {14,0,239}, {15,0,239}, {16,0,238},
 {17,0,237}, {18,0,236}, {19,0,234}, {20,0,233}, {21,0,232}, {22,0,231}, {23,0,231}, {24,0,230}, {25,0,229},
 {26,0,228}, {27,0,226}, {28,0,225}, {29,0,224}, {30,0,223}, {31,0,223}, {33,0,222}, {33,0,221}, {35,0,220},
 {35,0,218}, {37,0,217}, {37,0,216}, {39,0,215}, {39,0,215}, {41,0,214}, {41,0,213}, {43,0,212}, {43,0,210},
 {45,0,209}, {45,0,208}, {47,0,207}, {47,0,207}, {49,0,206}, {49,0,205}, {51,0,204}, {51,0,202}, {53,0,201},
 {53,0,200}, {55,0,199}, {55,0,199}, {57,0,198}, {57,0,197}, {59,0,196}, {59,0,194}, {61,0,193}, {61,0,192},
 {63,0,191}, {63,0,191}, {64,0,190}, {66,0,189}, {67,0,188}, {67,0,186}, {68,0,185}, {70,0,184}, {71,0,183},
 {71,0,183}, {72,0,182}, {74,0,181}, {75,0,180}, {75,0,178}, {76,0,177}, {78,0,176}, {79,0,175}, {79,0,175},
 {80,0,174}, {82,0,173}, {83,0,172}, {83,0,170}, {84,0,169}, {86,0,168}, {87,0,167}, {87,0,167}, {88,0,166},
 {90,0,165}, {91,0,164}, {91,0,162}, {92,0,161}, {94,0,160}, {95,0,159}, {95,0,159}, {96,0,158}, {98,0,157},
 {99,0,156}, {99,0,154}, {100,0,153}, {102,0,152}, {103,0,151}, {103,0,151}, {104,0,150}, {106,0,149}, {107,0,148},
{107,0,146}, {108,0,145}, {110,0,144}, {111,0,143}, {111,0,143}, {112,0,142}, {114,0,141}, {115,0,140}, {115,0,138},
{116,0,137}, {118,0,136}, {119,0,135}, {119,0,135}, {120,0,134}, {122,0,133}, {123,0,132}, {123,0,130}, {124,0,129},
{126,0,128}, {127,0,127}, {127,0,127}, {128,0,126}, {129,0,124}, {130,0,123}, {132,0,123}, {133,0,122}, {134,0,120},
{135,0,119}, {135,0,119}, {136,0,118}, {137,0,116}, {138,0,115}, {140,0,115}, {141,0,114}, {142,0,112}, {143,0,111},
{143,0,111}, {144,0,110}, {145,0,108}, {146,0,107}, {148,0,107}, {149,0,106}, {150,0,104}, {151,0,103}, {151,0,103},
{152,0,102}, {153,0,100}, {154,0, 99}, {156,0, 99}, {157,0, 98}, {158,0, 96}, {159,0, 95}, {159,0, 95}, {160,0, 94},
{161,0, 92}, {162,0, 91}, {164,0, 91}, {165,0, 90}, {166,0, 88}, {167,0, 87}, {167,0, 87}, {168,0, 86}, {169,0, 84},
{170,0, 83}, {172,0, 83}, {173,0, 82}, {174,0, 80}, {175,0, 79}, {175,0, 79}, {176,0, 78}, {177,0, 76}, {178,0, 75},
{180,0, 75}, {181,0, 74}, {182,0, 72}, {183,0, 71}, {183,0, 71}, {184,0, 70}, {185,0, 68}, {186,0, 67}, {188,0, 67},
{189,0, 66}, {190,0, 64}, {191,0, 63}, {191,0, 63}, {192,0, 61}, {193,0, 61}, {194,0, 59}, {196,0, 59}, {197,0, 57},
{198,0, 57}, {199,0, 55}, {199,0, 55}, {200,0, 53}, {201,0, 53}, {202,0, 51}, {204,0, 51}, {205,0, 49}, {206,0, 49},
{207,0, 47}, {207,0, 47}, {208,0, 45}, {209,0, 45}, {210,0, 43}, {212,0, 43}, {213,0, 41}, {214,0, 41}, {215,0, 39},
{215,0, 39}, {216,0, 37}, {217,0, 37}, {218,0, 35}, {220,0, 35}, {221,0, 33}, {222,0, 33}, {223,0, 31}, {223,0, 30},
{224,0, 29}, {225,0, 28}, {226,0, 27}, {228,0, 26}, {229,0, 25}, {230,0, 24}, {231,0, 23}, {231,0, 22}, {232,0, 21},
{233,0, 20}, {234,0, 19}, {236,0, 18}, {237,0, 17}, {238,0, 16}, {239,0, 15}, {239,0, 14}, {240,0, 13}, {241,0, 12},
{242,0, 11}, {244,0, 10}, {245,0,  9}, {246,0,  8}, {247,0,  7}, {247,0,  6}, {248,0,  5}, {249,0,  4}, {250,0,  3},
{252,0,  2}, {253,0,  1}, {254,0,  0}, {255,0,  0}};

    int i = v * 256;
    if (i > 255) i = 255;
    if (i < 0) i = 0;
    *r = ramp[i][0];
    *g = ramp[i][1];
    *b = ramp[i][2];
}

void cmap_tipsy(double v, int *r, int *g, int *b)
{
    if ( v < 0.0 ) v = 0.0;
    else if ( v > 1.0 ) v = 1.0;

    v *= 6;

    if      ( v < 1 ) { *r=0;           *g=0;               *b=255*v;       }
    else if ( v < 2 ) { *r=255*(v-1.0); *g=0;               *b=255;         }
    else if ( v < 3 ) { *r=255;         *g=0;               *b=255*(3.0-v); }
    else if ( v < 4 ) { *r=255;         *g=255*(v-3.0)*0.5; *b=0;           }
    else if ( v < 5 ) { *r=255;         *g=255;             *b=255*(v-5.0); }
    else              { *r=255;         *g=255;             *b=255;         }
}

void cmap_grey(double v, int *r, int *g, int *b)
{
    if ( v < 0.0 ) v = 0.0;
    else if ( v > 1.0 ) v = 1.0;

    *r =
    *g =
    *b = 255*v;
}
