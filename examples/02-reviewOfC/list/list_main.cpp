#include <iostream>

typedef struct List {
	int item_;
	struct List *next_;
	List():item_(0), next_(0){}
} List;

void insert_list(List **l, int x)
{
     List *p = new List;
	p->item_ = x;
	p->next_ = *l;
	*l = p;
}

List *print_list(List *l)
{
	if (l == 0) 
		return(0);
	else
		std::cout << l->item_;
	return( print_list(l->next_) );
}


int main()
{
	List *l = new List;
	insert_list(&l, 2);
	//insert_list(&l, 4);

	print_list(l);
    return 0;
}