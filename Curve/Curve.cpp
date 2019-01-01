/*
 * Copyright (c) 2013 Clément Bœsch
 * VapourSynth port by HolyWu
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <string>

#include <VapourSynth.h>
#include <VSHelper.h>

struct CurveData {
    VSNodeRef * node;
    const VSVideoInfo * vi;
    bool process[3];
    std::unique_ptr<uint16_t[]> graph[4];
};

struct keypoint {
    double x, y;
    std::shared_ptr<keypoint> next;
};

static inline void parsePoints(const char * s, std::shared_ptr<keypoint> & points, const int scale) {
    char * p = const_cast<char *>(s);
    std::shared_ptr<keypoint> last;

    // construct a linked list based on the key points string
    while (p && *p) {
        std::shared_ptr<keypoint> point = std::make_shared<keypoint>();

        point->x = std::strtod(p, &p);
        if (p && *p)
            p++;

        point->y = std::strtod(p, &p);
        if (p && *p)
            p++;

        if (point->x < 0. || point->x > 1. || point->y < 0. || point->y > 1.)
            throw std::string{ "invalid key point coordinates, x and y must be in the [0;1] range" };

        if (!points)
            points = point;

        if (last) {
            if (static_cast<int>(last->x * scale + 0.5) >= static_cast<int>(point->x * scale + 0.5))
                throw std::string{ "key point coordinates are too close from each other or not strictly increasing on the x-axis" };
            last->next = point;
        }

        last = point;
    }

    if (points && !points->next)
        throw std::string{ "only one point is defined, this is unlikely to behave as you expect" };
}

static inline int getNumberOfPoints(const keypoint * d) noexcept {
    int n = 0;
    while (d) {
        n++;
        d = d->next.get();
    }
    return n;
}

/**
 * Natural cubic spline interpolation
 * Finding curves using Cubic Splines notes by Steven Rauch and John Stockie
 */
static inline void interpolate(const keypoint * points, uint16_t * VS_RESTRICT y, const int lutSize, const int scale) {
    const keypoint * point = points;
    double xPrev = 0.;

    const int n = getNumberOfPoints(points); // number of splines

    if (n == 0) {
        for (int i = 0; i < lutSize; i++)
            y[i] = i;
        return;
    }

    double (*matrix)[3] = reinterpret_cast<double (*)[3]>(calloc(n, sizeof(*matrix)));
    double * h = reinterpret_cast<double *>(malloc((n - 1) * sizeof(*h)));
    double * r = reinterpret_cast<double *>(calloc(n, sizeof(*r)));
    if (!matrix || !h || !r)
        throw std::string{ "malloc failure (matrix/h/r)" };

    // h(i) = x(i+1) - x(i)
    int i = -1;
    for (point = points; point; point = point->next.get()) {
        if (i != -1)
            h[i] = point->x - xPrev;
        xPrev = point->x;
        i++;
    }

    // right-side of the polynomials, will be modified to contain the solution
    point = points;
    for (i = 1; i < n - 1; i++) {
        const double yp = point->y;
        const double yc = point->next->y;
        const double yn = point->next->next->y;
        r[i] = 6 * ((yn - yc) / h[i] - (yc - yp) / h[i - 1]);
        point = point->next.get();
    }

    constexpr int BD = 0; // sub  diagonal (below main)
    constexpr int MD = 1; // main diagonal (center)
    constexpr int AD = 2; // sup  diagonal (above main)

    // left side of the polynomials into a tridiagonal matrix
    matrix[0][MD] = matrix[n - 1][MD] = 1;
    for (i = 1; i < n - 1; i++) {
        matrix[i][BD] = h[i - 1];
        matrix[i][MD] = 2 * (h[i - 1] + h[i]);
        matrix[i][AD] = h[i];
    }

    // tridiagonal solving of the linear system
    for (i = 1; i < n; i++) {
        const double den = matrix[i][MD] - matrix[i][BD] * matrix[i - 1][AD];
        const double k = den ? 1. / den : 1.;
        matrix[i][AD] *= k;
        r[i] = (r[i] - matrix[i][BD] * r[i - 1]) * k;
    }
    for (i = n - 2; i >= 0; i--)
        r[i] = r[i] - matrix[i][AD] * r[i + 1];

    point = points;

    // left padding
    for (i = 0; i < static_cast<int>(point->x * scale + 0.5); i++)
        y[i] = std::min(std::max(static_cast<int>(point->y * scale + 0.5), 0), scale);

    // compute the graph with x=[x0..xN]
    i = 0;
    while (point->next) {
        const double yc = point->y;
        const double yn = point->next->y;

        const double a = yc;
        const double b = (yn - yc) / h[i] - h[i] * r[i] / 2. - h[i] * (r[i + 1] - r[i]) / 6.;
        const double c = r[i] / 2.;
        const double d = (r[i + 1] - r[i]) / (6. * h[i]);

        const int xStart = static_cast<int>(point->x * scale + 0.5);
        const int xEnd = static_cast<int>(point->next->x * scale + 0.5);
        for (int x = xStart; x <= xEnd; x++) {
            const double xx = static_cast<double>(x - xStart) / scale;
            const double yy = a + b * xx + c * xx * xx + d * xx * xx * xx;
            y[x] = std::min(std::max(static_cast<int>(yy * scale + 0.5), 0), scale);
        }

        point = point->next.get();
        i++;
    }

    // right padding
    for (i = static_cast<int>(point->x * scale + 0.5); i < lutSize; i++)
        y[i] = std::min(std::max(static_cast<int>(point->y * scale + 0.5), 0), scale);

    free(matrix);
    free(h);
    free(r);
}

template<typename T>
static void filter(const VSFrameRef * src, VSFrameRef * dst, const CurveData * const VS_RESTRICT d, const VSAPI * vsapi) noexcept {
    for (int plane = 0; plane < d->vi->format->numPlanes; plane++) {
        if (d->process[plane]) {
            const int width = vsapi->getFrameWidth(src, plane);
            const int height = vsapi->getFrameHeight(src, plane);
            const int stride = vsapi->getStride(src, plane) / sizeof(T);
            const T * srcp = reinterpret_cast<const T *>(vsapi->getReadPtr(src, plane));
            T * VS_RESTRICT dstp = reinterpret_cast<T *>(vsapi->getWritePtr(dst, plane));

            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++)
                    dstp[x] = d->graph[plane][srcp[x]];

                srcp += stride;
                dstp += stride;
            }
        }
    }
}

static void VS_CC curveInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi) {
    CurveData * d = static_cast<CurveData *>(*instanceData);
    vsapi->setVideoInfo(d->vi, 1, node);
}

static const VSFrameRef *VS_CC curveGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    const CurveData * d = static_cast<const CurveData *>(*instanceData);

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrameRef * src = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrameRef * fr[] = { d->process[0] ? nullptr : src, d->process[1] ? nullptr : src, d->process[2] ? nullptr : src };
        const int pl[] = { 0, 1, 2 };
        VSFrameRef * dst = vsapi->newVideoFrame2(d->vi->format, d->vi->width, d->vi->height, fr, pl, src, core);

        if (d->vi->format->bytesPerSample == 1)
            filter<uint8_t>(src, dst, d, vsapi);
        else
            filter<uint16_t>(src, dst, d, vsapi);

        vsapi->freeFrame(src);
        return dst;
    }

    return nullptr;
}

static void VS_CC curveFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    CurveData * d = static_cast<CurveData *>(instanceData);
    vsapi->freeNode(d->node);
    delete d;
}

static void VS_CC curveCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    std::unique_ptr<CurveData> d = std::make_unique<CurveData>();
    int err;

    d->node = vsapi->propGetNode(in, "clip", 0, nullptr);
    d->vi = vsapi->getVideoInfo(d->node);

    try {
        if (!isConstantFormat(d->vi) ||
            (d->vi->format->sampleType == stInteger && d->vi->format->bitsPerSample > 16) ||
            d->vi->format->sampleType == stFloat)
            throw std::string{ "only constant format 8-16 bit integer input supported" };

        const int preset = int64ToIntS(vsapi->propGetInt(in, "preset", 0, &err));

        const int numCurve = vsapi->propNumElements(in, "curve");

        const int m = vsapi->propNumElements(in, "planes");

        for (int i = 0; i < 3; i++)
            d->process[i] = (m <= 0);

        for (int i = 0; i < m; i++) {
            const int n = int64ToIntS(vsapi->propGetInt(in, "planes", i, nullptr));

            if (n < 0 || n >= d->vi->format->numPlanes)
                throw std::string{ "plane index out of range" };

            if (d->process[n])
                throw std::string{ "plane specified twice" };

            d->process[n] = true;
        }

        if (preset < 0 || preset > 10)
            throw std::string{ "preset must be 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, or 10" };

        if (numCurve > d->vi->format->numPlanes)
            throw std::string{ "more curves given than there are planes" };

        const char * curve[4] = {};

        for (int i = 0; i < numCurve; i++)
            curve[i] = vsapi->propGetData(in, "curve", i, nullptr);

        curve[3] = vsapi->propGetData(in, "master", 0, &err);

        if (preset == 1) {
            if (!(curve[0] && *curve[0]))
                curve[0] = "0.129/1 0.466/0.498 0.725/0";
            if (!(curve[1] && *curve[1]))
                curve[1] = "0.109/1 0.301/0.498 0.517/0";
            if (!(curve[2] && *curve[2]))
                curve[2] = "0.098/1 0.235/0.498 0.423/0";
        } else if (preset == 2) {
            if (!(curve[0] && *curve[0]))
                curve[0] = "0/0 0.25/0.156 0.501/0.501 0.686/0.745 1/1";
            if (!(curve[1] && *curve[1]))
                curve[1] = "0/0 0.25/0.188 0.38/0.501 0.745/0.815 1/0.815";
            if (!(curve[2] && *curve[2]))
                curve[2] = "0/0 0.231/0.094 0.709/0.874 1/1";
        } else if (preset == 3) {
            if (!(curve[3] && *curve[3]))
                curve[3] = "0/0 0.5/0.4 1/1";
        } else if (preset == 4) {
            if (!(curve[3] && *curve[3]))
                curve[3] = "0/0 0.149/0.066 0.831/0.905 0.905/0.98 1/1";
        } else if (preset == 5) {
            if (!(curve[3] && *curve[3]))
                curve[3] = "0/0 0.4/0.5 1/1";
        } else if (preset == 6) {
            if (!(curve[3] && *curve[3]))
                curve[3] = "0/0 0.305/0.286 0.694/0.713 1/1";
        } else if (preset == 7) {
            if (!(curve[3] && *curve[3]))
                curve[3] = "0/0 0.286/0.219 0.639/0.643 1/1";
        } else if (preset == 8) {
            if (!(curve[3] && *curve[3]))
                curve[3] = "0/1 1/0";
        } else if (preset == 9) {
            if (!(curve[3] && *curve[3]))
                curve[3] = "0/0 0.301/0.196 0.592/0.6 0.686/0.737 1/1";
        } else if (preset == 10) {
            if (!(curve[0] && *curve[0]))
                curve[0] = "0/0.11 0.42/0.51 1/0.95";
            if (!(curve[1] && *curve[1]))
                curve[1] = "0/0 0.50/0.48 1/1";
            if (!(curve[2] && *curve[2]))
                curve[2] = "0/0.22 0.49/0.44 1/0.8";
        }

        const int lutSize = 1 << d->vi->format->bitsPerSample;
        const int scale = lutSize - 1;
        std::shared_ptr<keypoint> points[4];

        for (int i = 0; i < 4; i++) {
            d->graph[i] = std::make_unique<uint16_t[]>(lutSize);
            parsePoints(curve[i], points[i], scale);
            interpolate(points[i].get(), d->graph[i].get(), lutSize, scale);
        }

        if (curve[3]) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < lutSize; j++)
                    d->graph[i][j] = d->graph[3][d->graph[i][j]];
            }
        }
    } catch (const std::string & error) {
        vsapi->setError(out, ("Curve: " + error).c_str());
        vsapi->freeNode(d->node);
        return;
    }

    vsapi->createFilter(in, out, "Curve", curveInit, curveGetFrame, curveFree, fmParallel, 0, d.release(), core);
}

//////////////////////////////////////////
// Init

VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin) {
    configFunc("com.holywu.curve", "curve", "Apply color adjustments using curves", VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("Curve",
                 "clip:clip;"
                 "preset:int:opt;"
                 "curve:data[]:opt;"
                 "master:data:opt;"
                 "planes:int[]:opt;",
                 curveCreate, nullptr, plugin);
}
