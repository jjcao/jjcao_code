#ifndef DGAL_TYPE_MANIP_H
#define DGAL_TYPE_MANIP_H

#include <DGAL/config.h>

DGAL_BEGIN_NAMESPACE

////////////////////////////////////////////////////////////////////////////////
// class template Int2Type
// Converts each integral constant into a unique type
// Invocation: Int2Type<v> where v is a compile-time constant integral
// Defines 'value', an enum that evaluates to v
////////////////////////////////////////////////////////////////////////////////

    template <int v>
    struct Int2Type
    {
        enum { value = v };
    };

DGAL_END_NAMESPACE

#endif //DGAL_TYPE_MANIP_H
