//
//  File:       CBLuts.cpp
//
//  Function:   Utilities for colour-blind modelling and LUT construction
//
//  Copyright:  Andrew Willmott 2018
//

#include "CBLuts.h"

#include <math.h>
#include <assert.h>

using namespace CBLut;

// --- Colour-blind support ---------------------------------------------------

// LMS colour space, models human eye response: https://en.wikipedia.org/wiki/LMS_color_space

// https://ixora.io/projects/colorblindness/color-blindness-simulation-research/
// More recent version of original approach below, uses more up-to-date LMS
// conversion, different approach to colour constraints, and observes Tritanope
// conversion has an issue in that it appears to have been derived by ensuring
// blue remains constant rather than red or green...

const Mat3f CBLut::kLMSFromRGB =
{
    0.31399022f,    0.63951294f,    0.04649755f,
    0.15537241f,    0.75789446f,    0.08670142f,
    0.01775239f,    0.10944209f,    0.87256922f,
};

const Mat3f CBLut::kRGBFromLMS =
{
    5.47221206f,   -4.64196010f,    0.16963708f,
    -1.1252419f,    2.29317094f,   -0.16789520f,
    0.02980165f,   -0.19318073f,    1.16364789f,
};

const Mat3f CBLut::kLMSProtanope =      /// Protanope: red sensitivity is greatly reduced, reds/yellows appear darker (1% men).
{
    0, 1.05118294f, -0.05116099f,
    0, 1, 0,
    0, 0, 1,
};

const Mat3f CBLut::kLMSDeuteranope =    /// Deuteranope: green sensivitity is greatly reduced, no brightness issues (1% men)
{
    1, 0, 0,
    0.9513092f, 0,  0.04866992f,
    0, 0, 1,
};

const Mat3f CBLut::kLMSTritanope =      /// Tritanope: blue sensitivity greatly reduced (0.003% population)
{
    1, 0, 0,
    0, 1, 0,
    -0.86744736f, 1.86727089f, 0
};

namespace
{
    inline float& elt(      Vec3f& v, int i) { return (&v.x)[i]; } 
//  inline Vec3f& row(      Mat3f& m, int i) { return (&m.x)[i]; } 
    inline float  elt(const Vec3f& v, int i) { return (&v.x)[i]; } 
    inline Vec3f  row(const Mat3f& m, int i) { return (&m.x)[i]; } 
    inline Vec3f  col(const Mat3f& m, int i) { return { (&m.x.x)[i], (&m.y.x)[i], (&m.z.x)[i] }; } 

    inline Vec3f operator+(Vec3f a, Vec3f b) { return { a.x + b.x, a.y + b.y, a.z + b.z}; }
    inline Vec3f operator-(Vec3f a, Vec3f b) { return { a.x - b.x, a.y - b.y, a.z - b.z}; }
    inline Vec3f operator*(float s, Vec3f a) { return { s   * a.x, s   * a.y, s   * a.z}; }
    
    inline float dot      (Vec3f a, Vec3f b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
    inline Vec3f pow      (Vec3f v, float p) { return { powf(v.x, p), powf(v.y, p), powf(v.z, p) }; }

    inline Vec3f operator*(const Mat3f& m, const Vec3f& v) { return Vec3f { dot(m.x, v), dot(m.y, v), dot(m.z, v) }; }
}

// "Digital Video Colourmaps for Checking the Legibility of Displays by Dichromats", Viénot et al.
//
// Example implementations: http://www.daltonize.org,
// Unfortunately now dead (sublinks redirect to a new site), copy(?) here:
//   https://github.com/joergdietrich/daltonize/blob/master/daltonize.py

namespace
{
    // Note that unlike kLMSFromRGB, LMS are weighted, e.g., red -> (17.8, 3.4, 0.02), blue ->  
    const Mat3f kLMSFromRGBV =
    {
         { 17.8824f,   43.5161f,   4.11935f },
         { 3.45565f,   27.1554f,   3.86714f },
         { 0.0299566f,  0.184309f, 1.46709f },
    };

    const Mat3f kRGBFromLMSV =
    {
        {  0.080944447900f, -0.13050440900f,  0.116721066f },
        { -0.010248533500f,  0.05401932660f, -0.113614708f },
        { -0.000365296938f, -0.00412161469f,  0.693511405f },
    };

    // These transforms to LMS colours simulate particular forms of colour blindness
    const Mat3f kLMSProtanopeV =      /// Protanope: red sensitivity is greatly reduced, reds/yellows appear darker (1% men).
    {
         { 0.0f, 2.02344f, -2.52581f, },
         { 0.0f, 1.0f,      0.0f,     },
         { 0.0f, 0.0f,      1.0f      },
    };

    const Mat3f kLMSDeuteranopeV =    /// Deuteranope: green sensivitity is greatly reduced, no brightness issues (1% men)
    {
         { 1.0f,      0.0f, 0.0f,      },
         { 0.494207f, 0.0f, 1.24827f,  },
         { 0.0f,      0.0f, 1.0f       },
    };

    const Mat3f kLMSTritanopeV =      /// Tritanope: blue sensitivity greatly reduced (0.003% population)
    {
         1.0f,       0.0f,      0.0f,
         0.0f,       1.0f,      0.0f,
         -0.395913f, 0.801109f, 0.0f
    };
}

namespace
{
    // From Onur Fidaner, Poliang Lin, and Nevran Ozguven. http://scien.stanford.edu/class/psych221/ projects/05/ofidaner/project report.pdf.
    // Unfortunately SCIEN have seen fit to break the links on their website, and there seems to be no author copy.
    // Update: a copy can be found here for now: https://github.com/joergdietrich/daltonize/blob/master/doc/project_report.pdf
    // The matrix values can be found here: https://github.com/joergdietrich/daltonize/blob/master/doc/conv_img.m.
    // Their precise values aren't discussed or justified in the paper.
    const Mat3f kDaltonErrorToDeltaP =
    {
        { 0.0f, 0.0f, 0.0f, },
        { 0.7f, 1.0f, 0.0f, },
        { 0.7f, 0.0f, 1.0f  },
    };
    const Mat3f kDaltonErrorToDeltaD =
    {
        { 1.0f, 0.7f, 0.0f, },
        { 0.0f, 0.0f, 0.0f, },
        { 0.0f, 0.7f, 1.0f  },
    };
    const Mat3f kDaltonErrorToDeltaT =
    {
        { 1.0f, 0.0f, 0.7f, },
        { 0.0f, 1.0f, 0.7f, },
        { 0.0f, 0.0f, 0.0f  },
    };

    Vec3f SimulateV(Vec3f rgb, const Mat3f& lmsTransform)
    {
        Vec3f lms = lmsTransform * (kLMSFromRGBV * rgb);
        return kRGBFromLMSV * lms;
    }
}

Vec3f CBLut::Daltonise(Vec3f rgb, tLMS lmsType, float strength)
{
    // Daltonisation: take delta from original RGB
    // + use to shift colors towards visible spectrum
    Vec3f rgbSim;
    Vec3f rgbDelta;
    
    switch (lmsType)
    {
    case kL:
        rgbSim = SimulateV(rgb, kLMSProtanopeV);
        rgbDelta = kDaltonErrorToDeltaP * (strength * (rgb - rgbSim));
        break;
    case kM:
        rgbSim = SimulateV(rgb, kLMSDeuteranopeV);
        rgbDelta = kDaltonErrorToDeltaD * (strength * (rgb - rgbSim));
        break;
    case kS:
        rgbSim = SimulateV(rgb, kLMSTritanopeV);
        rgbDelta = kDaltonErrorToDeltaT * (strength * (rgb - rgbSim));
        break;
    }
    
    return rgb + rgbDelta;
}

namespace
{
    const Mat3f kLMSSimulate =      // P/D/T simulation amalgamated into one matrix
    {
         0,            1.05118294f, -0.05116099f,
         0.9513092f,   0,            0.04866992f,
        -0.86744736f,  1.86727089f,  0
    };

    const Mat3f kNCDeltaBrighten =     // vanilla transfer error to remaining channels, like Daltonise approach but in lms space. amount=2.5 gets closish match
    {
        0, 1, 1,
        1, 0, 1,
        1, 1, 0,
    };

    const Mat3f kNCDeltaContrast =     // vanilla increase contrast between remaining channels. Not very controllable.
    {
         0, +1, +1,
        -1,  0, -1,
        +1, -1,  0,
    };

    const Mat3f kNCDeltaMix1 =          // mix of contrast/brighten
    {
        0.0f, 0.2f, 0.8f,
        0.2f, 0.0f, 0.2f,
        0.8f, 0.8f, 0.0f,
    };

    // const Mat3f kLMSDaltonError = kLMSFromRGB * kDaltonErrorToDeltaP * kRGBFromLMS;
    const Mat3f kNCDalton =         // The Fidaner RGB delta transformed to LMS space.
    {
        { 1.90957534f, -0.771573185f, 0.0281965993f },
        { 2.38503432f, -1.023173690f, 0.0739354342f },
        { 3.66449642f, -3.108514070f, 1.1135983500f },
    };

    const Mat3f kNCDeltaRecip =         // trans(1/kLMSSimulate)
    {
        0,             1.05118299f, -1.15280771f,
        0.951309144f,  0,            0.535540938f,
        -19.5461426f,  20.5465717f,  0,
    };

    const Mat3f kNCDeltaRecipAbs =      // abs(trans(1/kLMSSimulate))
    {
         0,              1.05118299f,  1.15280771f,
         0.951309144f,   0,            0.535540938f,
        19.5461426f,    20.5465717f,   0,
    };

    const Mat3f kNCDeltaInv =     // inv(kSimulate)
    {
          0.672f, 0.706f, -0.378f,
          0.312f, 0.328f,  0.360f,
        -13.133f, 6.741f,  7.393f
    };
}

Vec3f CBLut::Simulate(Vec3f rgb, tLMS lmsType, float strength)
{
    Vec3f lms = kLMSFromRGB * rgb;

    float&      eltx   = elt(lms, lmsType);                     // affected channel
    const float simElt = dot(row(kLMSSimulate, lmsType), lms);  // 'sim' weighted combo of the other two.

    eltx += strength * (simElt - eltx);

    rgb = kRGBFromLMS * lms;

    return rgb;
}

Vec3f CBLut::Correct(Vec3f rgb, tLMS lmsType, float strength)
{
    const Vec3f lms = kLMSFromRGB * rgb;

    const float orgElt = elt(lms, lmsType);                     // original value of affected channel
    const float simElt = dot(row(kLMSSimulate, lmsType), lms);  // simulated full-strength value
    const float error  = strength * (orgElt - simElt);          // error

    float mc = strength * strength; // How much to use strategy 1: redistributing error into other channels in a way that shifts hue
    float ms = 1.0f - strength;     // How much to use stragegy 2: simply brighten affected channel

    Vec3f amount3Recip = { -0.25f, -0.3f, -0.07f };    // tuning values for redistribution
    float amount = elt(amount3Recip, lmsType);

    Vec3f correct = mc * amount * col(kNCDeltaRecip, lmsType);
    elt(correct, lmsType) = ms * 2.0f;

    Vec3f lmsCorrect = lms + error * correct;
    
    rgb = kRGBFromLMS * lmsCorrect;

    return rgb;
}


// --- RGB LUT support ---------------------------------------------------------

namespace
{
    constexpr float kGamma = 2.2f;

    inline uint8_t ToU8(float f)
    {
        if (f <= 0.0f)
            return 0;
        if (f >= 1.0f)
            return 255;

        return uint8_t(f * 255.0f + 0.5f);
    }

    inline uint8_t ToU8u(float f)   // 0-256 variant used for LUT construction
    {
        if (f <= 0.0f)
            return 0;
        if (f >= 1.0f)
            return 255;

        return uint8_t(f * 256.0f);
    }
}

RGBA32 CBLut::ToRGBA32(Vec3f c)
{
    c = pow(c, 1.0f / kGamma);
    RGBA32 result;

    result.c[0] = ToU8(c.x);
    result.c[1] = ToU8(c.y);
    result.c[2] = ToU8(c.z);
    result.c[3] = 255;

    return result;
}

RGBA32 CBLut::ToRGBA32u(Vec3f c)
{
    c = pow(c, 1.0f / kGamma);
    RGBA32 result;

    result.c[0] = ToU8u(c.x);
    result.c[1] = ToU8u(c.y);
    result.c[2] = ToU8u(c.z);
    result.c[3] = 255;

    return result;
}

Vec3f CBLut::FromRGBA32(RGBA32 rgb)
{
    Vec3f c = { rgb.c[0] / 255.0f, rgb.c[1] / 255.0f, rgb.c[2] / 255.0f };
    return pow(c, kGamma);
}

Vec3f CBLut::FromRGBA32u(RGBA32 rgb)
{
    Vec3f c = { rgb.c[0] / 256.0f, rgb.c[1] / 256.0f, rgb.c[2] / 256.0f };
    return pow(c, kGamma);
}
    

void CBLut::CreateIdentityLUT(RGBA32 rgbLUT[kLUTSize][kLUTSize][kLUTSize])
{
    constexpr int scale  = 256 / kLUTSize;
    constexpr int offset = scale / 2;

    for (int i = 0; i < kLUTSize; i++)
    for (int j = 0; j < kLUTSize; j++)
    for (int k = 0; k < kLUTSize; k++)
    {
        RGBA32& p = rgbLUT[i][j][k];

        p.c[0] = k * scale + offset;
        p.c[1] = j * scale + offset;
        p.c[2] = i * scale + offset;
        p.c[3] = 255;
    }
}

#define EXTRAPOLATE_LUT 1

void CBLut::ApplyLUT(RGBA32 rgbLUT[kLUTSize][kLUTSize][kLUTSize], int n, const RGBA32 dataIn[], RGBA32 dataOut[])
{
    constexpr int lutShift = kLUTBits;
    constexpr int lutSize  = 1 << lutShift;
    constexpr int fShift   = 8 - lutShift;
    constexpr int fHalf    = 1 << (fShift - 1);
    constexpr int fMask    = (1 << fShift) - 1;

    for (int i = 0; i < n; i++)
    {
        const uint8_t* ci = dataIn[i].c;

        int co[3] = { ci[0] + fHalf,   ci[1] + fHalf,   ci[2] + fHalf   };
        int i1[3] = { co[0] >> fShift, co[1] >> fShift, co[2] >> fShift };
        int i0[3] = { i1[0] - 1,       i1[1] - 1,       i1[2] - 1       };
        int s [3] = { co[0] & fMask,   co[1] & fMask,   co[2] & fMask   };

        for (int j = 0; j < 3; j++)
        {
            if (i0[j] < 0)
            {
                i0[j]++;
            #ifdef EXTRAPOLATE_LUT
                i1[j]++;
                s [j] -= 8;
            #endif
            }
            else
            if (i1[j] >= lutSize)
            {
                i1[j]--;
            #ifdef EXTRAPOLATE_LUT
                i0[j]--;
                s [j] += 8;
            #endif
            }

            assert(0 <= i0[j] && i0[j] < kLUTSize);
            assert(0 <= i1[j] && i1[j] < kLUTSize);
        }
        
        RGBA32 lutC0 = rgbLUT[i0[2]][i0[1]][i0[0]];
        RGBA32 lutC1 = rgbLUT[i1[2]][i1[1]][i1[0]];

        int ch0 = (((8 - s[0]) * lutC0.c[0] + s[0] * lutC1.c[0])) >> 3;
        int ch1 = (((8 - s[1]) * lutC0.c[1] + s[1] * lutC1.c[1])) >> 3;
        int ch2 = (((8 - s[2]) * lutC0.c[2] + s[2] * lutC1.c[2])) >> 3;

    #ifdef EXTRAPOLATE_LUT
        ch0 = ch0 < 0 ? 0 : ch0 > 255 ? 255 : ch0;
        ch1 = ch1 < 0 ? 0 : ch1 > 255 ? 255 : ch1;
        ch2 = ch2 < 0 ? 0 : ch2 > 255 ? 255 : ch2;
    #endif

        assert(0 <= ch0 && ch0 <= 255);
        assert(0 <= ch1 && ch1 <= 255);
        assert(0 <= ch2 && ch2 <= 255);

        dataOut[i].c[0] = ch0;
        dataOut[i].c[1] = ch1;
        dataOut[i].c[2] = ch2;
        dataOut[i].c[3] = 255;
    }
}

void CBLut::ApplyLUTNoLerp(RGBA32 rgbLUT[kLUTSize][kLUTSize][kLUTSize], int n, const RGBA32 dataIn[], RGBA32 dataOut[])
{
    constexpr int fShift = 8 - kLUTBits;

    for (int i = 0; i < n; i++)
    {
        const uint8_t* ci = dataIn[i].c;

        dataOut[i] = rgbLUT[ci[2] >> fShift][ci[1] >> fShift][ci[0] >> fShift];
    }
}

// --- Mono LUT support --------------------------------------------------------

void CBLut::ApplyMonoLUT(const RGBA32 monoLUT[256], int n, const RGBA32 dataIn[], RGBA32 dataOut[], int channel)
{
    if (channel < 0)
    {
        for (int i = 0; i < n; i++)
        {
            Vec3f c = FromRGBA32(dataIn[i]);    // now linear
            float lumD65 = dot(Vec3f{0.2126f, 0.7152f, 0.0722f}, c);
            
            uint8_t lumU8 = ToU8(pow(lumD65, 1.0f / kGamma));    // lookup tables are in gamma space
        
            dataOut[i] = monoLUT[lumU8];
        }
        return;
    }
    
    for (int i = 0; i < n; i++)
        dataOut[i] = monoLUT[dataIn[i].c[channel]];
}
    
