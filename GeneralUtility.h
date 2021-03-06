#ifndef __GENERALUTILITY_H_CMC__
#define __GENERALUTILITY_H_CMC__

// Use iprint
#ifndef iprint
#define iprint(name) myprinter(#name, (name))
#endif
template <typename TType>
inline void myprinter (string name, const TType& value)
{
    cout << name << endl;
    cout << value << endl;
}

template <typename TypeName>
inline void mycheck (const TypeName& condition, string message="")
{
    if (!bool(condition))
    {
        cout << message << endl;
        throw;
    }
}

#endif
