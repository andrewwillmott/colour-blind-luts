//
//  File:       CBLutGen.cpp
//
//  Function:   Utilities for colour-blind modelling and LUT construction
//
//  Copyright:  Andrew Willmott 2018
//

#define _CRT_SECURE_NO_WARNINGS

#include "CBLuts.h"

#include "ColourMaps.h"

#include "stb_image_mini.h"

#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#ifdef _MSC_VER
    #define strlcpy(d, s, ds) strcpy_s(d, ds, s)
#endif

using namespace CBLut;

namespace
{
    inline Vec3f operator+(Vec3f a, Vec3f b) { return { a.x + b.x, a.y + b.y, a.z + b.z}; }
    inline Vec3f operator-(Vec3f a, Vec3f b) { return { a.x - b.x, a.y - b.y, a.z - b.z}; }
    
    inline Vec3f& operator*=(Vec3f& v, float s) { v.x *= s; v.y *= s; v.z *= s; return v; }
    
    inline float dot      (Vec3f a, Vec3f b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
    inline Vec3f operator*(Vec3f a, Vec3f b) { return { a.x * b.x, a.y * b.y, a.z * b.z}; }

    inline Vec3f operator * (const Mat3f& m, const Vec3f& v)
    {
        return Vec3f
        {
            dot(m.x, v),
            dot(m.y, v),
            dot(m.z, v),
        };
    }

    inline Vec3f ClampUnit(Vec3f c)
    {
        return { 
        c.x < 0.0f ? 0.0f : c.x > 1.0f ? 1.0f : c.x, 
        c.y < 0.0f ? 0.0f : c.y > 1.0f ? 1.0f : c.y, 
        c.z < 0.0f ? 0.0f : c.z > 1.0f ? 1.0f : c.z };
    }

    Vec3f LMSError(Vec3f rgb, const Mat3f& lmsTransform)
    {
        Vec3f lms = (kLMSFromRGB * rgb);
        Vec3f lmsSim = lmsTransform * lms;
        return lmsSim - lms;
    }

    Vec3f RGBError(Vec3f c, tLMS lmsType, float strength)
    {
        Vec3f sc = Simulate(c, lmsType, strength); 
        return c - sc;
    }
    
    Vec3f LMSSwap(Vec3f rgb, tLMS ch)
    {
        Vec3f lms = kLMSFromRGB * rgb;
        float t;

        switch (ch)
        {
        case kL:
            t = lms.x; lms.x = lms.y; lms.y = t; 
            break;
        case kM:
            t = lms.y; lms.y = lms.z; lms.z = t; 
            break;
        case kS:
            t = lms.z; lms.z = lms.x; lms.x = t; 
            break;
        };
        
        rgb = kRGBFromLMS * lms;
        return rgb;
    }

    Vec3f RemapLToS(Vec3f rgb)
    {
        Vec3f lmsP = kLMSFromRGB * rgb;
        Vec3f lmsSimP = kLMSProtanope * lmsP;
        
        float error = lmsP.x - lmsSimP.x; 

        Vec3f lmsS = lmsSimP;
        lmsS.z += 10 * error;
        
        return kRGBFromLMS * lmsS;
    }

    Vec3f RemapMToS(Vec3f rgb)
    {
        Vec3f lmsM = kLMSFromRGB * rgb;
        Vec3f lmsSimM = kLMSDeuteranope * lmsM;
        
        float error = lmsM.x - lmsSimM.x; 

        Vec3f lmsS = lmsSimM;
        lmsS.z += 10 * error;
        
        return kRGBFromLMS * lmsS;
    }

    template<class T> void CreateLUT(T xform, RGBA32 rgbLUT[kLUTSize][kLUTSize][kLUTSize])
    {
        constexpr int scale  = 256 / kLUTSize;
        constexpr int offset = scale / 2;

        for (int i = 0; i < kLUTSize; i++)
        for (int j = 0; j < kLUTSize; j++)
        for (int k = 0; k < kLUTSize; k++)
        {
            // Vec3f c{ (k + 0.5f) / kLUTSize, (j + 0.5f) / kLUTSize, (i + 0.5f) / kLUTSize };
            RGBA32 identity = { uint8_t(k * scale + offset), uint8_t(j * scale + offset), uint8_t(i * scale + offset), 255 };

            Vec3f c = FromRGBA32u(identity);

            c = xform(c);

            rgbLUT[i][j][k] = ToRGBA32u(c);
        }
    }

    template<class T> void Transform(T xform, int n, const RGBA32 dataIn[], RGBA32 dataOut[])
    {
        for (int i = 0; i < n; i++)
        {
            Vec3f c = FromRGBA32(dataIn[i]);

            c = xform(c);

            dataOut[i] = ToRGBA32(c);
        }
    }

    template<class T> inline void PerformOp(T xform, RGBA32 rgbLUT[kLUTSize][kLUTSize][kLUTSize], int n, const RGBA32 dataIn[], RGBA32 dataOut[])
    {
        if (dataOut)
            Transform(xform, n, dataIn, dataOut);
        else
            CreateLUT(xform, rgbLUT);
    }
}


namespace
{
    // Main rgb LUT processing routines
    enum tCBType
    {
        kIdentity,

        kProtanope,
        kDeuteranope,
        kTritanope,

        kAll,
    };

    enum tImageOp 
    {
        kSimulate,
        kError,
        kDaltonise,
        kCorrect,
        kDaltoniseSimulate,
        kCorrectSimulate,
        kPassThrough,
    };

    void CreateImage(tImageOp op, tCBType cbType, float strength, int w, int h, const RGBA32* dataIn, const char* dataInName, bool noLUT)
    {
        if (cbType == kAll)
        {
            CreateImage(op, kProtanope,   strength, w, h, dataIn, dataInName, noLUT);
            CreateImage(op, kDeuteranope, strength, w, h, dataIn, dataInName, noLUT);
            CreateImage(op, kTritanope,   strength, w, h, dataIn, dataInName, noLUT);
            return;
        };

        tLMS lmsType = kL;
        char filename[256] = "";

        if (dataIn)
            snprintf(filename, sizeof(filename), "%s_", dataInName);

        switch (cbType)
        {
        case kIdentity:
            strcat(filename, "identity");
            break;
        case kProtanope:
            strcat(filename, "protanope");
            lmsType = kL;
            break;
        case kDeuteranope:
            strcat(filename, "deuteranope");
            lmsType = kM;
            break;
        case kTritanope:
            strcat(filename, "tritanope");
            lmsType = kS;
            break;
        default:
            return;
        }

        RGBA32 rgbaLUT[kLUTSize][kLUTSize][kLUTSize];
        RGBA32* dataOut = 0;
        int n = w * h;
        
        if (noLUT && dataIn) 
            dataOut = new RGBA32[n];
        
        switch (op)
        {
        case kSimulate:
            PerformOp([lmsType, strength](Vec3f c){ return Simulate(c, lmsType, strength); }, rgbaLUT, n, dataIn, dataOut);
            strcat(filename, "_simulate");
            break;
        case kError:
            PerformOp([lmsType, strength](Vec3f c){ return RGBError(c, lmsType, strength); }, rgbaLUT, n, dataIn, dataOut);
            strcat(filename, "_error");
            break;
        case kDaltonise:
            PerformOp([lmsType, strength](Vec3f c) { return Daltonise(c, lmsType, strength); }, rgbaLUT, n, dataIn, dataOut);
            strcat(filename, "_daltonise");
            break;
        case kCorrect:
            PerformOp([lmsType, strength](Vec3f c) { return Correct(c, lmsType, strength); }, rgbaLUT, n, dataIn, dataOut);
            strcat(filename, "_correct");
            break;
        case kDaltoniseSimulate:
            PerformOp([lmsType, strength](Vec3f c) { return Simulate(ClampUnit(Daltonise(c, lmsType, strength)), lmsType, strength); }, rgbaLUT, n, dataIn, dataOut);
            strcat(filename, "_simulate_daltonised");
            break;
        case kCorrectSimulate:
            PerformOp([lmsType, strength](Vec3f c) { return Simulate(ClampUnit(Correct(c, lmsType, strength)), lmsType, strength); }, rgbaLUT, n, dataIn, dataOut);
            strcat(filename, "_simulate_corrected");
            break;
        case kPassThrough:
            if (dataOut)
                PerformOp([](Vec3f c) { return c; }, rgbaLUT, n, dataIn, dataOut);
            else
                CreateIdentityLUT(rgbaLUT);
            break;
        };

        if (dataIn && !dataOut)
        {
            dataOut = new RGBA32[n];

            ApplyLUT(rgbaLUT, n, dataIn, dataOut);
        }

        if (dataOut)
        {
            strcat(filename, ".png");
            printf("Saving %s\n", filename);
            stbi_write_png(filename, w, h, 4, dataOut, 0);

            delete[] dataOut;
        }
        else
        {
            strcat(filename, "_lut.png");
            printf("Saving %s\n", filename);
            stbi_write_png(filename, kLUTSize * kLUTSize, kLUTSize, 4, rgbaLUT, 0);
        }
    }

    void CreateImage(const RGBA32* rgbaLUT, int w, int h, const RGBA32* dataIn)
    {
        int n = w * h;
        RGBA32* dataOut = new RGBA32[n];

        ApplyLUT(* (RGBA32 (*)[kLUTSize][kLUTSize][kLUTSize]) (RGBA32*) rgbaLUT, w * h, dataIn, dataOut);
        
        char filename[256] = "apply_lut";
        
        printf("Saving %s\n", filename);
        strcat(filename, ".png");

        stbi_write_png(filename, w, h, 4, dataOut, 0);

        delete[] dataOut;
    }
}


namespace
{
    // Mono LUT processing
    struct cMonoLUTEntry { const char* name; const uint8_t (*lut)[4]; } kMonoLUTs[]  = 
    {
        "cividis", kCividisLUT,
        "viridis", kViridisLUT,
        "magma",   kMagmaLUT  ,
        "inferno", kInfernoLUT,
        "plasma",  kPlasmaLUT ,
    };

    void PrintMonoLUT(const char* name, const RGBA32 monoLUT[256])
    {
        printf("const unsigned char k%s[256][4] =\n{\n", name);

        for (int i = 0; i < 256; i++)
        {
            const uint8_t* c = monoLUT[i].c;
            printf("    %3d, %3d, %3d, %3d,\n", c[0], c[1], c[2], c[3]);
        } 

        printf("};\n");
    }

    void CreateImageWithMonoLUT(const RGBA32 monoLUT[256], const char* lutName, int w, int h, const RGBA32* dataIn, const char* dataName, int channel)
    {
        RGBA32* dataOut = 0;

        if (dataIn)
        {
            dataOut = new RGBA32[w * h];
            ApplyMonoLUT(monoLUT, w * h, dataIn, dataOut, channel);
        }
        else
        {
            w = 256;
            h = 8;
            dataOut = new RGBA32[w * h];
            for (int i = 0; i < h; i++)
                memcpy(dataOut + i * w, monoLUT, w * sizeof(RGBA32));
        }
        
        char filename[256];
        
        if (dataIn)
            snprintf(filename, sizeof(filename), "%s_%s.png", dataName, lutName);
        else
            snprintf(filename, sizeof(filename), "%s_lut.png", lutName);
        
        printf("Saving %s\n", filename);
        stbi_write_png(filename, w, h, 4, dataOut, 0);

        delete[] dataOut;
    }
}

namespace
{
    int Help(const char* command)
    {
        printf
        (
            "%s <options> <operations>\n"
            "\n"
            "Options:\n"
            "  -h        : this help\n"
            "  -f <path> : set image to process rather than emitting lut\n"
            "  -p        : emit protanope image or lut\n"
            "  -d        : emit deuteranope image or lut\n"
            "  -t        : emit tritanope image or lut\n"
            "  -a        : emit image or lut for all the above types (default)\n"
            "  -m <str>  : specify strength of colour blindness to correct for. Default = 1 (affected channel is completely lost.)\n" 
            "  -n        : directly transform input image rather than using a LUT\n"
            "  -g[LMS]   : swap LM/MS/LS channels of input image before processing\n"
            "  -r[LM]    : remap L or M channels to S, converting a prot/deuter test image to tritanope.\n"
            "\n"
            "Operations:\n"
            "  -s        : simulate given type of colour-blindness\n"
            "  -x        : daltonise (Fidaner) for given type of colour-blindness\n"
            "  -X        : daltonise for and then simulate given type of colour-blindness\n"
            "  -y        : correct for given type of colour-blindness\n"
            "  -Y        : correct for and then simulate given type of colour-blindness\n"
            "  -e        : error between original colour and simulated version\n"
            "  -i        : emit identity image or lut (for testing)\n"
            "  -l <path> : apply the given LUT to source (requires -f)\n"
            "\n"
            "  -c <name> [<channel>] : apply given greyscale lut: cividis, viridis (cb-savvy). magma, inferno, plasma (standard)\n"
            "                          'name' can also be the path of a 256-wide LUT in image form\n"
            "                          if channel is supplied, it is used to index the lut, otherwise sRGB/D65 luminance is used\n"
            "\nExample:\n"
            "  %s -f image.png -p -sxy\n"
            "      # emit simulated, daltonised, and corrected version of image.png for protanopia only.\n"
            , command, command
        );

        return 0;
    }

    void GetFileName(char* buffer, size_t bufferSize, const char* path)
    {
        const char* lastSlash = strrchr(path, '/');
        if (!lastSlash)
            lastSlash = strrchr(path, '\\');

        if (lastSlash)
            strlcpy(buffer, lastSlash + 1, bufferSize);
        else
            strlcpy(buffer, path, bufferSize);

        char* lastDot = strrchr(buffer, '.');

        if (lastDot)
            *lastDot = 0;
    }
}

int main(int argc, const char* argv[])
{
    const char* command = argv[0];
    argv++; argc--;

    if (argc == 0)
        return Help(command);

    tCBType cbType = kAll;
    int w;
    int h;
    RGBA32* dataIn = 0;
    char dataInName[256] = "unknown";
    float strength = 1.0f;
    bool noLUT = false;

    // Options
    while (argc > 0 && argv[0][0] == '-')
    {
        const char* option = argv[0] + 1;
        argv++; argc--;

        while (option[0])
        {
            switch (option[0])
            {
            case 'h':
            case '?':
                return Help(command);

            case 'c':
                {
                    const RGBA32* lutTable = 0;
                    const char*   lutName  = 0;
                    int           channel  = -1;

                    if (!(argc > 0 && argv[0][0] != '-'))
                    {
                        fprintf(stderr, "Expecting lut argument\n");
                        return -1;
                    }

                    for (cMonoLUTEntry& entry : kMonoLUTs)
                        if (strcmp(argv[0], entry.name) == 0)
                        {
                            lutTable = (const RGBA32*) entry.lut;
                            lutName  = entry.name;
                            break;
                        }

                    if (!lutTable)
                    {
                        int lw, lh;
                        lutTable = (RGBA32*) stbi_load(argv[0], &lw, &lh, 0, 4);

                        static char lutNameStore[256];
                        GetFileName(lutNameStore, sizeof(lutNameStore), argv[0]);
                        lutName = lutNameStore;

                        if (!lutTable)
                        {
                            fprintf(stderr, "Unknown mono LUT or file not found: %s\n", argv[0]);
                            return -1;
                        }

                        if (lw != 256)
                        {
                            fprintf(stderr, "Expecting mono LUT width of 256\n");
                            return -1;
                        }
                    }

                    argv++; argc--;

                    if (argc > 0 && argv[0][0] != '-')
                    {
                        channel = atoi(argv[0]);
                        argv++; argc--;
                    }
                    
                    CreateImageWithMonoLUT(lutTable, lutName, w, h, dataIn, dataInName, channel);
                        // PrintMonoLUT(lutName, lutTable);
                }
                break;

            case 'f':
                if (argc <= 0)
                    return fprintf(stderr, "Expecting filename with -f\n");

                dataIn = (RGBA32*) stbi_load(argv[0], &w, &h, 0, 4);
                
                if (!dataIn)
                {
                    fprintf(stderr, "Couldn't read %s\n", argv[0]);
                    return -1;
                }

                GetFileName(dataInName, sizeof(dataInName), argv[0]);

                argv++; argc--;
                break;

            case 'F':
                {
                    // Create a swatch that varies horizontally only in L, for
                    // protanope correction testing.
                    w = 256;
                    h = 256;
                    dataIn = new RGBA32[w * h];
                    strcpy(dataInName, "swatch");

                    RGBA32* p = dataIn;
                    for (int y = 0; y < h; y++)
                    for (int x = 0; x < w; x++)
                    {
                        Vec3f lms;

                    #if 0
                        // tiles
                        lms.x = ((x & ~0x1F) + 16) / float(w);
                        lms.y = ((y & ~0x1F) + 16) / float(h);
                    #else
                        lms.x = (x + 0.5f) / float(w);
                        lms.y = (y + 0.5f) / float(h);
                    #endif

                        lms.z = 1.0f - lms.y;

                        // in LMS space, L and M are usually close to the same (because of their large overlap),
                        // and slight differences lead to red or green. Thus to stay within RGB gamut we must
                        // heavily restrict their range. S on the other hand is quite independent.
                        Vec3f remapLMS = Vec3f{ 0.46f, 0.45f, 0.25f  } + Vec3f{ 0.08f, 0.1f, 0.5f } * lms;
                        remapLMS *= 0.75;
                        Vec3f rgb = kRGBFromLMS * remapLMS;
                        
                        assert(rgb.x >= 0.0f && rgb.x <= 1.0f);
                        assert(rgb.y >= 0.0f && rgb.y <= 1.0f);
                        assert(rgb.z >= 0.0f && rgb.z <= 1.0f);
                        
                        RGBA32 c = ToRGBA32(rgb);

                        (*p++) = c;
                    }
                }
                break;

            case 'm':
                if (argc <= 0)
                    return fprintf(stderr, "Expecting strength for -m <float>\n");
                strength = (float) atof(argv[0]);
                argv++; argc--;
                break;

            case 's':
                CreateImage(kSimulate,          cbType, strength, w, h, dataIn, dataInName, noLUT);
                break;

            case 'e':
                CreateImage(kError,             cbType, strength, w, h, dataIn, dataInName, noLUT);
                break;

            case 'x':
                CreateImage(kDaltonise,         cbType, strength, w, h, dataIn, dataInName, noLUT);
                break;
            case 'X':
                CreateImage(kDaltoniseSimulate, cbType, strength, w, h, dataIn, dataInName, noLUT);
                break;

            case 'y':
                CreateImage(kCorrect,           cbType, strength, w, h, dataIn, dataInName, noLUT);
                break;
            case 'Y':
                CreateImage(kCorrectSimulate,   cbType, strength, w, h, dataIn, dataInName, noLUT);
                break;

            case 'i':
                CreateImage(kPassThrough, kIdentity, strength, w, h, dataIn, dataInName, noLUT);
                break;

            case 'g':
                if (option[1] == 'l' or option[1] == 'L')
                    Transform([](Vec3f c){ return LMSSwap(c, kL); }, w * h, dataIn, dataIn);
                else if (option[1] == 'm' or option[1] == 'M')
                    Transform([](Vec3f c){ return LMSSwap(c, kM); }, w * h, dataIn, dataIn);
                else
                    Transform([](Vec3f c){ return LMSSwap(c, kS); }, w * h, dataIn, dataIn);
                option++;

            case 'r':
                if (option[1] == 'm' or option[1] == 'M')
                    Transform([](Vec3f c){ return RemapMToS(c); }, w * h, dataIn, dataIn);
                else
                    Transform([](Vec3f c){ return RemapLToS(c); }, w * h, dataIn, dataIn);
                option++;
                break;

            case 'p':
                cbType = kProtanope;
                break;

            case 'd':
                cbType = kDeuteranope;
                break;

            case 't':
                cbType = kTritanope;
                break;

            case 'a':
                cbType = kAll;
                break;
                
            case 'n':
                noLUT = true;
                break;

            case 'l':
                if (argc <= 0)
                    return fprintf(stderr, "Expecting filename with -l\n");

                if (!dataIn)
                    return fprintf(stderr, "No input file to apply lut to\n");

                int lw, lh;
                RGBA32* lut = (RGBA32*) stbi_load(argv[0], &lw, &lh, 0, 4);
                
                if (!lut)
                {
                    fprintf(stderr, "Couldn't read RGB LUT %s\n", argv[0]);
                    return -1;
                }

                if (lw != kLUTSize * kLUTSize)
                {
                    fprintf(stderr, "Expecting RGB LUT width of %d\n", kLUTSize * kLUTSize);
                    return -1;
                }
                
                if (lh != kLUTSize)
                {
                    fprintf(stderr, "Expecting RGB LUT height of %d\n", kLUTSize);
                    return -1;
                }

                CreateImage(lut, w, h, dataIn);
                
                argv++; argc--;
                break;
            }
            option++;
        }
    }

    if (dataIn)
        stbi_image_free(dataIn);

    if (argc > 0)
    {
        fprintf(stderr, "Unrecognised arguments starting with %s\n", argv[0]);
        return -1;
    }
        
    return 0;
}
