#ifndef LINKED_LIST_H
#define LINKED_LIST_H

#include <src/box.h>

using namespace std;

template<class T> class linkedList
{
friend class box;
//protected:
public:
    T item;
    linkedList* next;
public:
    linkedList():
        next(0)
    {
//        cout << "using linked list default constructor" << endl;
    } // default constructor

    linkedList(T& item, linkedList* N=0):
        item(item),
        next(N)
    {
//        cout << "using linked list regular constructor" << endl;
    } // constructor

    linkedList(const linkedList& list): // copy constructor
        item(list.readItem()),
        next(list.next ? new linkedList(*list.next) : 0)
    {
//        cout << "using linked list copy constructor" << endl;
    } // copy constructor

    ~linkedList()
    {
        delete next;
        next = 0;

//        cout << "using linked list destructor" << endl;
    } // destructor

    const T& readItem() const // previously operator()()
    {
        return item;
    } // read item field

    const linkedList* readNext() const
    {
        return next;
    } // read next field

//    linkedList& last()
//    {
//        return next ? next->last() : *this;
//    } // last item

    int length() const
    {
        return next ? next->length() + 1 : 1;
    } // number of items

//    void append(T& t)
//    {
//        //cout << "&t = " << &t << endl;
//        last().next = new linkedList(t);
//    } // append a new item

    void insertNextItem(T& t)
    {
        next = new linkedList(t, next);
    } // inserts item in the second place

    void insertFirstItem(T& t)
    {
        next = new linkedList(item, next);
        item = t;
    } // insert a new item at the beginning ("push_front")

    const linkedList& operator=(const linkedList&);
    void dropNextItem();
    void dropFirstItem();
};

template<class T>
const linkedList<T>&
linkedList<T>::operator=(const linkedList<T>&L)
{
    if (this != &L)
    {
        item = L.readItem();
        if (next)
        {
            if (L.next)
            {
                *next = *L.next;
            }
            else
            {
                delete next;
                next = 0;
            }
        }
        else
            if (L.next) next = new linkedList(*L.next);
    }
    return *this;
} // assignment operator

template<class T>
void linkedList<T>::dropNextItem()
{
    if (next)
    {
        if (next->next)
        {
            linkedList<T>* keep = next;
            next = next->next;
            keep->item.~T();
        }
        else
        {
            delete next;
            next = 0;
        }
    }
    else
        cout << "error: cannot drop nonexisting next item" << endl;
} // drop the second item from the linked list

template<class T>
void linkedList<T>::dropFirstItem()
{
    if (next)
    {
        item = next->item;
        dropNextItem();
    }
    else
        cout << "error: cannot drop the only item" << endl;
} // drop the first item in the linked list



#endif // LINKED_LIST_H
