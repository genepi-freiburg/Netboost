#ifndef __RCPtr_hpp
#define __RCPtr_hpp

namespace RCPtrNamespace {

  struct tArrayFlagOn {
  };

  struct tArrayFlagOff {
  };

  template <typename T, typename tArrayFlag>
  class PerObj;

  template <typename T, typename tArrayFlag>
  class RCPtrBase {
    friend class RCPtrBase<T const, tArrayFlag>;
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef T& reference;
    typedef void iterator_category;
    typedef void difference_type;

  protected:
    RCPtrBase(T* p, bool responsible = true) :
      perObjPtr(new PerObj<T, tArrayFlag>(p, responsible)) {
    }
    
    // CCtor
    RCPtrBase(RCPtrBase const & other) throw() :
      perObjPtr(other.perObjPtr) {
      perObjPtr->inc();
    }

    ~RCPtrBase() throw() {
      if(perObjPtr->dec() == 0)
	delete perObjPtr;
    }

    RCPtrBase& operator=(RCPtrBase const & other) throw() {
      if(perObjPtr != other.perObjPtr) {
	if(perObjPtr->dec() == 0)
	  delete perObjPtr;
	perObjPtr = other.perObjPtr;
	perObjPtr->inc();
      }
      return *this;
    }

  public:
    T& operator*() const throw() {
      return perObjPtr->ref();
    }

    T* operator->() const throw() {
      return perObjPtr->ptr();
    }

    operator bool() const {
      return perObjPtr->ptr();
    }

    bool isSoleOwner() const {
      return perObjPtr->isSoleOwner();
    }

  private:
    PerObj<T, tArrayFlag>* perObjPtr;
  };


  template <typename T, typename tArrayFlag>
  class RCPtrBase<T const, tArrayFlag> {
    friend class RCPtrBase<T, tArrayFlag>;
  public:
    typedef T const value_type;
    typedef T const * pointer;
    typedef T const & reference;
    typedef void iterator_category;
    typedef void difference_type;

  protected:
    RCPtrBase(T const * p, bool responsible = true) :
      perObjPtr(new PerObj<T, tArrayFlag>(const_cast<T *>(p), responsible)) {
    }

    // CCtor
    RCPtrBase(RCPtrBase const & other) throw() :
      perObjPtr(other.perObjPtr) {
      perObjPtr->inc();
    }

    RCPtrBase(RCPtrBase<T, tArrayFlag> const & other) throw() :
      perObjPtr(other.perObjPtr) {
      perObjPtr->inc();
    }

    ~RCPtrBase() throw() {
      if(perObjPtr->dec() == 0)
	delete perObjPtr;
    }

  public:
    RCPtrBase& operator=(RCPtrBase const & other) throw() {
      if(perObjPtr != other.perObjPtr) {
	if(perObjPtr->dec() == 0)
	  delete perObjPtr;
	perObjPtr = other.perObjPtr;
	perObjPtr->inc();
      }
      return *this;
    }

    RCPtrBase& operator=(RCPtrBase<T, tArrayFlag> const & other) throw() {
      if(perObjPtr != other.perObjPtr) {
	if(perObjPtr->dec() == 0)
	  delete perObjPtr;
	perObjPtr = other.perObjPtr;
	perObjPtr->inc();
      }
      return *this;
    }

    T const & operator*() const throw() {
      return perObjPtr->ref();
    }

    T const * operator->() const throw() {
      return perObjPtr->ptr();
    }

    operator bool() const {
      return perObjPtr->ptr();
    }

    bool isSoleOwner() const {
      return perObjPtr->isSoleOwner();
    }

  private:
    PerObj<T, tArrayFlag>* perObjPtr;
  };


  template <typename T1, typename T2, typename A1, typename A2>
  bool operator==(RCPtrBase<T1, A1> const & l, RCPtrBase<T2, A2> const & r) {
    return (l.operator->() == r.operator->());
  }

  template <typename T1, typename T2, typename A>
  bool operator==(RCPtrBase<T1, A> const & l, T2* r) {
    return (l.operator->() == r);
  }

  template <typename T1, typename T2, typename A>
  bool operator==(T1* l, RCPtrBase<T1, A> const & r) {
    return (l == r.operator->());
  }

  template <typename T1, typename T2, typename A1, typename A2>
  bool operator!=(RCPtrBase<T1, A1> const & l, RCPtrBase<T2, A2> const & r) {
    return (l.operator->() != r.operator->());
  }

  template <typename T1, typename T2, typename A>
  bool operator!=(RCPtrBase<T1, A> const & l, T2* r) {
    return (l.operator->() != r);
  }

  template <typename T1, typename T2, typename A>
  bool operator!=(T1* l, RCPtrBase<T1, A> const & r) {
    return (l != r.operator->());
  }

  


  template <typename T, typename tArrayFlag>
  class PerObj {
  public:
    PerObj(T* p, bool responsible) throw() :
      p(p),
      count(1),
      responsible(responsible) {
    }

    ~PerObj() throw() {
      if (responsible)
	del(tArrayFlag());
    }
    
    void inc() throw() {
      ++count;
    }
    
    unsigned int dec() throw() {
      return --count;
    }
    
    T& ref() {
      return *p;
    }
    
    T* ptr() {
      return p;
    }

    bool isSoleOwner() const {
      return (count == 1);
    }
  private:
    void del(tArrayFlagOn) {
      delete[] p;
    }

    void del(tArrayFlagOff) {
      delete p;
    }
    
    
    PerObj(PerObj const&);
    void operator=(PerObj const &);
    
    T* p;
    unsigned int count;
    bool responsible;
  };


  struct DefaultToNull {
  };

  struct DefaultToNewObject {
  };
  
  template <typename T, typename DefaultPolicy = DefaultToNull, typename tArrayFlag = tArrayFlagOff>
  class RCPtr : public RCPtrBase<T, tArrayFlag> {
    typedef RCPtrBase<T, tArrayFlag> Base;
  public:
    // Ctor - initialize only with address of object created with new.
    //        do not initialize more than one RCPtr with the same address.
    RCPtr() :
      Base(GetDefault(DefaultPolicy()), true) {
    }

    RCPtr(T* p, bool responsible = true) :
      Base(p, responsible) {
    }

    RCPtr(RCPtr const & other) :
      Base(other) {
    }

    template <typename OT, typename tOArrayFlag>
    RCPtr(RCPtrBase<OT, tOArrayFlag> const & other) :
      Base(other) {
    }

    template <typename OT, typename tOArrayFlag>
    RCPtr& operator=(RCPtrBase<OT, tOArrayFlag> const  & other) {
      return static_cast<RCPtr&>(Base::operator=(other));
    }

    template <typename OT>
    RCPtr& operator=(OT* ptr) {
      return operator=(static_cast<RCPtr>(ptr));
    }

    ~RCPtr() throw() {
    }

    operator T*() const {
      return RCPtrBase<T, tArrayFlag>::operator->();
    }

  private:
    static T* GetDefault(DefaultToNull) {
      return 0;
    }

    static T* GetDefault(DefaultToNewObject) {
      return new T;
    }
  };
}

using RCPtrNamespace::RCPtr;


#endif
