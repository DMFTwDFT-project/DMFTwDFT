// @Copyright 2007 Kristjan Haule
// 

#include "CXX/Objects.hxx"
#include "CXX/Extensions.hxx"
#include <assert.h>
#include <algorithm>
#include <vector>
#include <deque>
#include <list>

class gutils_module : public Py::ExtensionModule<gutils_module>
{
public:

  gutils_module() : Py::ExtensionModule<gutils_module>( "gutils" )
  {
    add_varargs_method("compres", &gutils_module::ex_compres, "list[list[int]] compres(list[list[int]]) compresses lists");
    initialize( "This modeule provides utils written in C++ : ....." );
  }
  
  virtual ~gutils_module(){}

private:
  
  Py::Object ex_compres_orig( const Py::Tuple &a);
  Py::Object ex_compres( const Py::Tuple &a);
};

// required function initializes connection with Python
extern "C" void initgutils()
{ static gutils_module* gutils = new gutils_module;}

// symbol required for the debug version
extern "C" void initgutils_d()
{ initgutils(); }


void debug_check_ref_queue()
{
#ifdef Py_TRACE_REFS

    // create an element to find the queue
    Py::Int list_element;

    PyObject *p_slow = list_element.ptr();
    PyObject *p_fast = p_slow;

    do
    {
        assert( p_slow->_ob_next->_ob_prev == p_slow );
        assert( p_slow->_ob_prev->_ob_next == p_slow );


        p_slow = p_slow->_ob_next;
        p_fast = p_slow->_ob_next->_ob_next;

        assert( p_slow != p_fast );    
    }
    while( p_slow != list_element.ptr() );

#endif
}

static Py::List Union(const Py::List& data1, const Py::List& data2)
{
  // Takes a union of two lists
  Py::List res = data1;
  for (Py::List::const_iterator d=data2.begin(); d!=data2.end(); d++){
    bool d_in_res = false;
    for (Py::List::iterator dr=res.begin(); dr!=res.end(); dr++){
      if (*d == *dr){
	d_in_res=true;
	break;
      }
    }
    if (!d_in_res) res.append(*d);
  }
  return res;
}

bool Overlap(const Py::List& data1, const Py::List& data2)
{
  // Checks if there is any overlap between data1 and data2"
  for (Py::List::const_iterator d=data2.begin(); d!=data2.end(); d++){
    for (Py::List::const_iterator d1=data1.begin(); d1!=data1.end(); d1++){
      if (*d == *d1) return true;
    }
  }
  return false;
}

static Py::List Compres(Py::List& groups){
  bool loopsDone = true;
  while (loopsDone){
    loopsDone = false;
    for (int  i=0; i<groups.length(); i++){
      if (loopsDone) break;
      for (int j=i+1; j<groups.length(); j++){
	if (loopsDone) break;
	
	Py::List gi = static_cast<Py::List>(groups[i]);
	Py::List gj = static_cast<Py::List>(groups[j]);

	if (Overlap(gi,gj)){
	  groups[i] = Union(gi,gj);
	  
	  groups[j]=Py::List();
	  
	  loopsDone = true;
	}
      }
    }
  }
  Py::List ngroups;
  for (int ig=0; ig<groups.length(); ig++){
    Py::List gg = static_cast<Py::List>(groups[ig]);
    if (gg.length()>0){
      gg.sort();
      ngroups.append(gg);
    }
  }

  ngroups.sort();
    
  return ngroups;
}

Py::Object gutils_module::ex_compres_orig( const Py::Tuple &a){

  debug_check_ref_queue();
  
  Py::List lst(a[0]);
  return Compres(lst);
}

//List getSlice (int i, int j) const
//void setSlice (int i, int j, const Object& v)
//void append (const Object& ob)
//void insert (int i, const Object& ob)
//void sort ()
//void reverse ()
//void clear ()
//size_type size ()
//int length ()
//bool hasKey (const std::string& s) const
//bool hasKey (const Object& s) const

using namespace std;

template<class container>
bool Overlap(const container& data1, const container& data2)
{
  // Checks if there is any overlap between data1 and data2"
  for (typename container::const_iterator d=data2.begin(); d!=data2.end(); d++){
    for (typename container::const_iterator d1=data1.begin(); d1!=data1.end(); d1++){
      if (*d == *d1) return true;
    }
  }
  return false;
}

template<class container>
container Union(const container& data1, const container& data2)
{
  // Takes a union of two lists
  container result(data1);
  for (typename container::const_iterator d=data2.begin(); d!=data2.end(); d++){
    bool d_in_res = false;
    for (typename container::iterator dr=result.begin(); dr!=result.end(); dr++){
      if (*d == *dr){
	d_in_res=true;
	break;
      }
    }
    if (!d_in_res) result.push_back(*d);
  }
  return result;
}

// deque<int> Union(const deque<int>& data1, const deque<int>& data2)
// {
//   // Takes a union of two lists
//   deque<int> result(data1);
//   for (deque<int>::const_iterator d=data2.begin(); d!=data2.end(); d++){
//     bool d_in_res = false;
//     for (deque<int>::iterator dr=result.begin(); dr!=result.end(); dr++){
//       if (*d == *dr){
// 	d_in_res=true;
// 	break;
//       }
//     }
//     if (!d_in_res) result.push_back(*d);
//   }
//   return result;
// }

// bool Overlap(const deque<int>& data1, const deque<int>& data2)
// {
//   // Checks if there is any overlap between data1 and data2"
//   for (deque<int>::const_iterator d=data2.begin(); d!=data2.end(); d++){
//     bool d_in_d1 = false;
//     for (deque<int>::const_iterator d1=data1.begin(); d1!=data1.end(); d1++){
//       if (*d == *d1) return true;
//     }
//   }
//   return false;
// }

deque<deque<int> > Compres(vector<deque<int> >& groups){
  bool loopsDone = true;
  while (loopsDone){
    loopsDone = false;
    for (int i=0; i<groups.size(); i++){
      if (loopsDone) break;
      for (int j=i+1; j<groups.size(); j++){
	if (loopsDone) break;

	if (Overlap(groups[i],groups[j])){
	  groups[i] = Union(groups[i],groups[j]);
	  
	  groups[j].clear();
	  
	  loopsDone = true;
	}
      }
    }
  }
  deque<deque<int> > ngroups;
  for (int ig=0; ig<groups.size(); ig++){
    if (groups[ig].size()>0){
      sort(groups[ig].begin(),groups[ig].end());
      ngroups.push_back(groups[ig]);
    }
  }

  //ngroups.sort();
    
  return ngroups;
}


typedef list<list<int> > dbllist;

bool cmpr(list<int>& a, list<int>& b)
{
  return *(a.begin()) < *(b.begin());
}

dbllist Compress(dbllist& groups){
  bool loopsDone = true;
  dbllist::iterator i,j;
  while (loopsDone){
    loopsDone = false;
    for (i=groups.begin(); i!=groups.end(); ++i){
      if (loopsDone) break;
      j=i;
      for (++j; j!=groups.end(); ++j){
	if (loopsDone) break;

	if (Overlap(*i,*j)){
	  *i = Union(*i,*j);
	  
	  j->clear();
	  
	  loopsDone = true;
	}
      }
    }
  }
  dbllist ngroups;
  dbllist::iterator ig;
  for (ig=groups.begin(); ig!=groups.end(); ig++){
    if ((*ig).size()>0){
      (*ig).sort();
      ngroups.push_back(*ig);
    }
  }

  ngroups.sort(cmpr);
  
  return ngroups;
}

Py::Object gutils_module::ex_compres( const Py::Tuple &a){

  debug_check_ref_queue();
  
  Py::List lst(a[0]);

  
  dbllist inp;
  Py::List y;
  for(Py::List::size_type i=0; i < lst.length(); ++i){
    y = lst[i];
    list<int> yc;
    for(Py::List::size_type j=0; j < y.length(); ++j){
      int r = Py::Int(y[j]);
      yc.push_back(r);
    }
    inp.push_back(yc);
  }

  dbllist rs = Compress(inp);

  Py::List res;
  for (list<list<int> >::const_iterator i=rs.begin(); i!=rs.end(); i++){
    Py::List tres;
    for (list<int>::const_iterator j=(*i).begin(); j!=(*i).end(); j++){
      tres.append(Py::Int(*j));
    }
    res.append(tres);
  }
  return res;
}
