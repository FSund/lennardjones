#ifndef INLINES_H
#define INLINES_H

inline uint sub2ind3d(const uvec3 &subs, const uvec3 &subLimits)
{
    /* converts 3d subscripts (indexes) to a linear index */
    return subs(0) + subs(1)*subLimits(0) + subs(2)*subLimits(0)*subLimits(1);
}

inline uvec3 ind2sub3d(const uint &ind, const uvec3 &subLimits)
{
    /* converts linear index to 3d subscripts (indexes) */
    uvec3 sub;
    sub << ind%subLimits(0)
        << (ind/subLimits(0))%subLimits(1)
        << ((ind/subLimits(0))/subLimits(1))%subLimits(2);
    return sub;
}

//inline double lennard_jones(const vec3 &drvec)
//{
//    double dr2 = drvec(0)*drvec(0) + drvec(1)*drvec(1) + drvec(2)*drvec(2);
//    double dr6 = dr2*dr2*dr2;

//    return 24.0*(2.0 - dr6)/(dr6*dr6*dr2);
//}

#endif // INLINES_H
