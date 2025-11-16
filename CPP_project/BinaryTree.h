#ifndef BINARYTREE_H
#define BINARYTREE_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

template<typename T>
class BinaryTree{
    private:
        int _depth;
        std::vector<std::vector<T> > _tree;
    public:
        BinaryTree();
        void setDepth(int d);
        void setNode(int a, int b, T v);
        T getNode(int a, int b);
        void display();
};

// Template implementation must be in header file
template<typename T>
BinaryTree<T>::BinaryTree()
{
    _depth = 0;
}

template<typename T>
void BinaryTree<T>::setDepth(int d)
{
    _depth = d;
    _tree.resize(d+1);
    for(int i = 0 ; i<=d ; i++)
    {
        _tree[i].resize(i+1, T());  // Initialize with default value T()
    }
}

template<typename T>
void BinaryTree<T>::setNode(int a, int b, T v)
{
    _tree[a][b]=v;
}

template<typename T>
T BinaryTree<T>::getNode(int a, int b)
{
    return _tree[a][b];
}

template<typename T>
void BinaryTree<T>::display()
{
    for(int i = 0 ; i<=_depth; i++)
    {
        for(int j = 0 ; j<=i; j++)
        {
            std::cout << _tree[i][j] << "  ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n\n";
}

#endif