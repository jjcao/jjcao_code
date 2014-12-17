//===========================================================================
// GoTools - SINTEF Geometry Tools version 1.0.1
//
// GoTools module: CORE
//
// Copyright (C) 2000-2005 SINTEF ICT, Applied Mathematics, Norway.
//
// This program is free software; you can redistribute it and/or          
// modify it under the terms of the GNU General Public License            
// as published by the Free Software Foundation version 2 of the License. 
//
// This program is distributed in the hope that it will be useful,        
// but WITHOUT ANY WARRANTY; without even the implied warranty of         
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
// GNU General Public License for more details.                           
//
// You should have received a copy of the GNU General Public License      
// along with this program; if not, write to the Free Software            
// Foundation, Inc.,                                                      
// 59 Temple Place - Suite 330,                                           
// Boston, MA  02111-1307, USA.                                           
//
// Contact information: e-mail: tor.dokken@sintef.no                      
// SINTEF ICT, Department of Applied Mathematics,                         
// P.O. Box 124 Blindern,                                                 
// 0314 Oslo, Norway.                                                     
//
// Other licenses are also available for this software, notably licenses
// for:
// - Building commercial software.                                        
// - Building software whose source code you wish to keep private.        
//===========================================================================
#ifndef _GOFACTORY_H
#define _GOFACTORY_H

#include <map>
#include "GeomObject.h"
#include <boost/smart_ptr.hpp>


namespace Go
{
///\addtogroup geometry
///\{

    /// Abstract base class for Creators.  These are policy classes for generating
    /// objects derived from GeomObject.  This class is only intended for internal use 
    /// by the \ref Factory class.  The user should not have to worry about it.
    class Creator
    {
    public:
	virtual GeomObject* create() = 0;
	virtual ~Creator() {}
    };

    /// This is the concrete Creator class for generating GeomObject-derived classes.
    /// The class is only intended for internal use by the \ref Factory class.  The user
    /// should not have to worry about it.
    template <class T>
    class ConcreteCreator : public Creator
    {
    public:
	ConcreteCreator() {}
	GeomObject* create()
	{
	    return new T;
	}
    };
    
    /** This is the Factory for creating instances of GeomObjects of type
     *  requested by the user.  There will be one, global Factory object,
     *  called Go::global_factory_.  However, the user will not have to deal
     *  with this directly.  The user need only interact with the Factory 
     *  using the static member function 'createObject(ClassType)' and the
     *  template Register() function (or Registrator class, if the compiler
     *  misbehaves.
     */
    class Factory
    {
    public:
	~Factory()
	{
	    std::map<ClassType, Creator*>::iterator it;
	    for (it = themap_.begin(); it != themap_.end(); ++it)
		delete it->second;
	}

	/// Makes a new GeomObject instance of the specified ClassType and returns
	/// a pointer to it.  The user assumes ownership over the created object.
	/// \param class_type the class type of the object that the user wants to have
	///                   constructed.
	static GeomObject* createObject(ClassType class_type)
	{
	    return globalFactory()->doCreateObject(class_type);
	}

	/// Register a ClassType with the Factory.  This amounts to provide the Factory
	/// with the Creator object that is used to generate a new GeomObject of type
	/// ClassType.  Usually, the user would not want to call this function directly,
	/// but through the (template) function \ref Register() or by using a
	/// \ref Registrator
	/// \param class_type the ClassType for which we want the Factory to associate
	///                   a particular Creator.
	/// \param c pointer to the Creator object to be used by the Factory when creating
	///          objects of type ClassType.
	void registerClass(ClassType class_type, Creator* c)
	{
	    themap_[class_type] = c;
	}

	/// This function returns a pointer to the unique, global Factory.
	/// \return pointer to global Factory.
	static Factory* globalFactory()
	{
	    if (global_factory_.get() == 0)
		global_factory_ = boost::shared_ptr<Factory>(new Factory);
	    return global_factory_.get();
	}
    private:
	// private constructor.  User should not need to explicitly construct 
	// Factory.
	Factory()
	{
	}

	GeomObject* doCreateObject(ClassType class_type)
	{
	    std::map<ClassType, Creator*>::iterator it;
	    it = themap_.find(class_type);
	    if (it==themap_.end()) {
		THROW("Class type number " << class_type
		      << " is not registered.");
	    }
	    return it->second->create();
	}
	static boost::shared_ptr<Factory> global_factory_;
	std::map<ClassType, Creator*> themap_;
    };

    /// This function is used to register a class derived from GeomObject
    /// with the global Factory.  By using this function rather than
    /// \ref Factory::registerClass(), the user does not have to worry about the
    /// details of the Creator class.  
    /// To register a class \code DerivedClass \endcode, it should be 
    /// sufficient to run: \code Register<DerivedClass>() \endcode .
    template <class T>
    void Register()
    {
	Factory* f = Factory::globalFactory();
	ConcreteCreator<T>* c = new ConcreteCreator<T>;
	f->registerClass(T::classType(), c);
    }

    /** Missing doxygen documentation
     *
     */

    // @afr: I have no idea why, but VS6 refuses to do 

    /// On some compilators (ie., VS6), the \ref Register() function cannot
    /// be used directly.  To work around this, we wrap it in the constructor
    /// of a dummy class, which we call registrator.  With other words, in 
    /// order to register the class 'DerivedClass' (supposedly inheriting from
    /// GeomObject), we type: \code Registrator<DerivedClass> r \endcode.
    template <class T>
    class Registrator
    {
    public:
	Registrator()
	{
	    Factory* f = Factory::globalFactory();
	    ConcreteCreator<T>* c = new ConcreteCreator<T>;
	    f->registerClass(T::classType(), c);
	}
    };


///\}
} // namespace Go

#endif // _GOFACTORY_H


