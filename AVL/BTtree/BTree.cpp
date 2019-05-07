//增删查改
#include "iostream"
#include "string"

using namespace std;

typedef int ElementType;
typedef struct AVLNode
{
	ElementType data;
	struct AVLNode *pLeft;
	struct AVLNode *pRight;
	int Height;
}*Position, *AVLTree;

AVLTree MakeEmpty(AVLTree T)
{
	if (NULL != T)
	{
		MakeEmpty(T->pLeft);
		MakeEmpty(T->pRight);
		free(T);
	}
	return NULL;
}

Position Find(ElementType x, AVLTree T)
{
	if (NULL == T)
	{
		return NULL;
	}
	if (x < T->data)
	{
		return Find(x, T->pLeft);
	}
	else if (x > T->data)
	{
		return Find(x, T->pRight);
	}
	else
	{
		return T;
	}
}

Position FindMin(AVLTree T)
{
	if (NULL == T)
	{
		return NULL;
	}
	else if (NULL == T->pLeft)
	{
		return T;
	}
	else
	{
		return FindMin(T->pLeft);
	}
}

Position FindMax(AVLTree T)
{
	if (NULL != T)
	{
		while (NULL != T->pRight)
		{
			T = T->pRight;
		}
	}
	return T;
}

static int Height(Position p)
{
	if (NULL == p)
	{
		return -1;
	}
	else
	{
		return p->Height;
	}
}

static int Max(int h1, int h2)
{
	return h1 > h2 ? h1 : h2;
}
//LL右旋
static Position SingleRotateWithLeft(Position k2)//LL
{
	Position k1;
	k1 = k2->pLeft;
	k2->pLeft = k1->pRight;
	k1->pRight = k2;

	k2->Height = Max(Height(k2->pLeft), Height(k2->pRight)) + 1;
	k1->Height = Max(Height(k1->pLeft), Height(k1->pRight)) + 1;
	return k1;
}

//RR左旋
static Position SingleRotateWithRight(Position k1)
{
	Position k2;
	k2 = k1->pRight;
	k1->pRight = k2->pLeft;
	k2->pLeft = k1;
	k2->Height = Max(Height(k2->pLeft), Height(k2->pRight)) + 1;
	k1->Height = Max(Height(k1->pLeft), Height(k1->pRight)) + 1;
	return k2;
}

//LR B左旋 A右旋
static Position DoubleRotateLeft(Position k3)
{
	k3->pLeft = SingleRotateWithRight(k3->pLeft);
	return SingleRotateWithLeft(k3);
}

//RL B右旋 A左旋
static Position DoubleRotateRight(Position k4)
{
	k4->pRight = SingleRotateWithLeft(k4->pRight);
	return SingleRotateWithRight(k4);
}

AVLTree Insert(ElementType x, AVLTree T)
{
	if (NULL == T)
	{
		T = (AVLTree)malloc(sizeof(struct AVLNode));
		T->data = x;
		T->Height = 0;
		T->pLeft = T->pRight = NULL;
	}
	else if (x < T->data)
	{
		T->pLeft = Insert(x, T->pLeft);
		if ((Height(T->pLeft) - Height(T->pRight)) == 2)
		{
			if (x < T->pLeft->data)
			{
				T = SingleRotateWithLeft(T);
			}
			else
			{
				T = DoubleRotateLeft(T);
			}
		}
	}
	else if (x > T->data)
	{
		T->pRight = Insert(x, T->pRight);
		if ((Height(T->pRight)-Height(T->pLeft)) == 2)
		{
			if (x > T->pRight->data)
			{
				T = SingleRotateWithRight(T);
			}
			else
			{
				T = DoubleRotateRight(T);
			}
		}
	}
	T->Height = Max(Height(T->pLeft), Height(T->pRight)) + 1;
	return T;
}

AVLTree Rotate(AVLTree T)
{
	if (Height(T->pLeft) - Height(T->pRight) == 2)
	{
		if (Height(T->pLeft->pLeft) >= Height(T->pLeft->pRight))
			T = SingleRotateWithLeft(T);  // LL旋转  
		else
			T = DoubleRotateLeft(T);     // LR旋转  
	}
	if (Height(T->pRight) - Height(T->pLeft) == 2)
	{
		if (Height(T->pRight->pRight) >= Height(T->pRight->pLeft))
			T = SingleRotateWithRight(T);  // RR旋转  
		else
			T = DoubleRotateRight(T);     // RL旋转  
	}
	return T;
	
}

/*
* 首先定位要删除的节点，然后用该节点的右孩子的最左孩子替换该节点，
* 并重新调整以该节点为根的子树为AVL树，具体调整方法跟插入数据类似
* 删除处理在整体上耗费O(log n) 时间。
*/
AVLTree Delete(ElementType x, AVLTree T)
{
	if (T == NULL)
		return NULL;
	if (T->data == x)           // 要删除的 x 等于当前节点元素  
	{
		if (T->pRight == NULL)  // 若所要删除的节点 T 的右孩子为空,则直接删除  
		{
			AVLTree tmp = T;
			T = T->pLeft;
			free(tmp);
		}
		else                 /* 否则找到 T->right 的最左儿子代替 T */
		{
			AVLTree tmp = T->pRight;
			while (tmp->pLeft)
				tmp = tmp->pLeft;
			T->data = tmp->data;
			/* 对于替代后的T 即其字节点进行调整*/
			T->pRight = Delete(T->data, T->pRight);
			T->Height = Max(Height(T->pLeft), Height(T->pRight)) + 1;
		}
		return T;
	}
	else if (x > T->data)                       // 要删除的 x 大于当前节点元素，在T的右子树中查找删除  
	{
		T->pRight = Delete(x, T->pRight);
	}
	else                                       // 要删除的 x 小于当前节点元素，在T的左子树中查找删除  
	{
		T->pLeft = Delete(x, T->pLeft);
	}
	/*
	*   当删除元素后调整平衡
	*/
	T->Height = Max(Height(T->pLeft), Height(T->pRight)) + 1;
	if (T->pLeft != NULL)
		T->pLeft = Rotate(T->pLeft);
	if (T->pRight != NULL)
		T->pRight = Rotate(T->pRight);
	if (T)
		T = Rotate(T);
	return T;
}

ElementType Retrieve(Position P)
{
	return P->data;
}

/*
* 遍历输出
*/
void Display(AVLTree T)
{
	static int n = 0;
	if (NULL != T)
	{
		Display(T->pLeft);
		printf("[%d] ndata=%d \n", ++n, T->data);
		Display(T->pRight);
	}
}

#define N 15  
int main()
{
	AVLTree T = NULL;
	int i;
	int j = 0;
	T = MakeEmpty(NULL);
	for (i = 0; i < N; i++, j = (j + 7) % 50)
	{
		printf("j=%d \n", j);
		T = Insert(j, T);
	}
	puts("插入 4 \n");
	T = Insert(4, T);
	Display(T);
	for (i = 0; i < N; i += 2)
	{
		printf("delelte: %d \n", i);
		T = Delete(i, T);
	}
	printf("detele:\n");
	printf("height=%d \n", T->Height);
	Display(T);

	printf("Min is %d, Max is %d\n", Retrieve(FindMin(T)),
		Retrieve(FindMax(T)));
	getchar();
	return EXIT_SUCCESS;
	
}