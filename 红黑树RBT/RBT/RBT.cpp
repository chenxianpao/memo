// RBT.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include <stdio.h> 
#include <stdlib.h> 
#include <time.h>
#include "stdlib.h"
#include <string.h> 


/************************************************************************/
/* 
����������ʱ�临�Ӷ�ΪO��lgn�������в�����Ϊ��ά�ֺ�������������
1.ÿ�����Ҫô�Ǻ��Ҫô�Ǻڵġ�
2.������Ǻڵġ�
3.ÿ��Ҷ��㣨Ҷ��㼴ָ��β��NILָ���NULL��㣩���Ǻڵġ�
4.���һ������Ǻ�ģ���ô�����������Ӷ��Ǻڵġ�
5.������������ԣ��䵽Ҷ�����β��NILָ���ÿ��·����������ͬ��Ŀ�ĺڽ�㡣


���ֵ�ƽ����������Ч�ʣ��㷺����STL�С������Ǻ�������㷨ʵ�֣�һ��ѧϰ��
ƽʱ���Ҳ������STL��Ҳ�����о�����ô����Ԫ�صġ�
*/
/************************************************************************/
typedef int key_t;
typedef int data_t;

typedef enum color_t
{
	RED = 0,
	BLACK = 1
}color_t;

typedef struct rb_node_t
{
	struct rb_node_t *left, *right, *parent;
	key_t key;                  //����ǹؼ�ֵ
	data_t data;                //����Ǳ��
	color_t color;
}rb_node_t;

/*forward declaration */

/*������*/
rb_node_t *rb_insert(key_t key, data_t data, rb_node_t *root);
/*���ҽ��*/
rb_node_t *rb_search(key_t key, rb_node_t *root);
/*�ͷŽ��*/
rb_node_t *rb_erase(key_t key, rb_node_t *root);


static rb_node_t *rb_new_node(key_t key, data_t data)
{
	rb_node_t *node = (rb_node_t*)malloc(sizeof(struct rb_node_t));
	if (!node)
	{
		printf("malloc error!\\n");
		exit(-1);
	}
	node->key = key;
	node->data = data;
	return node;
}

/*-----------------------------------------------------------
|  node                  right

|  /\\         ==>         /\\
|  a right             node y

|      /\\              /\\

|      b y            a  b        
-----------------------------------------------------------*/
/* ����ת���� */
static rb_node_t *rb_rotate_left(rb_node_t *node, rb_node_t *root)
{

	rb_node_t *right = node->right;          //ָ��ָ��ָ��right

	if ((node->right = right->left))       //��right�����ӽ�� �ҽ� ��node���Һ���
	{
		right->left->parent = node;       //node��Ϊright���ӵĸ�ĸ���
	}

	right->left = node;               //node��Ϊright������

	if ((right->parent = node->parent))        //�ж�node��Ϊ�����
	{
		if (node == node->parent->right)       //nodeΪ���׽����Һ���     
		{
			node->parent->right = right;
		}
		else                    //nodeΪ���׽�������
		{
			node->parent->left = right;
		}
	}
	else
	{
		root = right;               //���nodeΪ�����    
	}
	node->parent = right;                    //right��Ϊnode�ĸ�ĸ
	return root;
}

/*-----------------------------------------------------------

|       node                 left
|        /\\                  /\\

|     left y  ==>          a  node
|      /\\                      /\\

|     a  b                    b  y //������������࣬�����Թ�

-----------------------------------------------------------*/
/* �������� */
static rb_node_t *rb_rotate_right(rb_node_t *node, rb_node_t *root)
{
	rb_node_t *left = node->left;            //ָ��ָ��ָ��left

	if ((node->left = left->right))           //��left���Һ��Ӹ�ֵ��node�����ӣ���left���Һ��Ӳ�Ϊ��
	{
		left->right->parent = node;       //node��Ϊleft�Һ��ӵĸ��׽��
	}
	left->right = node;              //node��Ϊleft���Һ���
	if ((left->parent = node->parent))         //��left�ĸ����ָ��node�ĸ���㣬���Ҵ�ʱ��node��Ϊ�����    
	{
		if (node == node->parent->right)       //nodeΪ�丸�����Һ���         
		{
			node->parent->right = left;
		}
		else                    //nodeΪ�丸��������
		{
			node->parent->left = left;
		}
	}
	else                        //nodeΪ�����  
	{
		root = left;
	}

	node->parent = left;             //left��Ϊnode�ĸ����
	return root;

}

//��������ҽ��
//-------------------------------------------------------------
//rb_search_auxiliary������
//rb_node_t * rb_search�������ҵ��Ľ��
//-------------------------------------------------------------
static rb_node_t *rb_search_auxiliary(key_t key, rb_node_t *root, rb_node_t **save)
{
	rb_node_t *node = root, *parent = NULL;
	int ret;
	while (node)
	{
		parent = node;
		ret = node->key - key;
		if (0 < ret)             //key С�ڵ�ǰ����key�����ԣ���������������
		{
			node = node->left;
		}
		else if (0 > ret)        //key ���ڵ�ǰ����key�����ԣ���������������
		{
			node = node->right;
		}
		else
		{
			return node;        //���ҳɹ�������
		}
	}
	if (save)                    //
	{
		*save = parent;
	}

	return NULL;
}

//��������rb_search_auxiliary���ҽṹ
rb_node_t *rb_search(key_t key, rb_node_t *root)
{
	return rb_search_auxiliary(key, root, NULL);
}

/* ������֮�����µ���������Ľṹ */

//����������������z��ʾ��ǰ��㣬p[z]��ʾ��ĸ��p[p[z]]��ʾ�游��y��ʾ����
static rb_node_t *rb_insert_rebalance(rb_node_t *node, rb_node_t *root)
{
	rb_node_t *parent, *gparent, *uncle, *tmp;
	//ѭ��������parentΪnode�ĸ�����Ҳ�Ϊ�գ����⣬node�ĸ�������ɫΪ��ɫʱ
	while ((parent = node->parent) && (parent->color == RED))
	{
		gparent = parent->parent;

		if (parent == gparent->left)                      //���游������Ϊ�����ʱ
		{                                  //��ʵ����������䣬�޷Ǿ�����˳���ӡ���ĸ���游�Ĺ�ϵ
			uncle = gparent->right;                      //��������ĸ������y���Ǹ�ĸ���Һ���

			if (uncle && uncle->color == RED)                     //���1��z��������y�Ǻ�ɫ�ģ�������Ϊ��ɫ�����뱣֤����y��ΪNULL��    
			{
				uncle->color = BLACK;                    //�Բ����£�a��b��c
				parent->color = BLACK;                   //a. ��node�ĸ�����������Ϳ��

				gparent->color = RED;                     //b. ���游���Ϳ��
				node = gparent;                      //c. �ѵ�ǰ���ָ���游��㣬��ָ��z����������������Ŷ

				//�������1ֻ������z�ĸ����Ϊ�游�����ӵ����
			}
			else                                   //���2��z��������y�Ǻ�ɫ��
			{
				if (parent->right == node)                    //��zΪ�Һ���
				{                           //�Բ����£�
					root = rb_rotate_left(parent, root);        //��ǰ���ĸ������Ϊ�µĵ�ǰ���
					tmp = parent;                   //���µĵ�ǰ���Ϊ֧������
					parent = node;
					node = tmp;
				}
				//���3��z������y�Ǻ�ɫ�ģ���ʱz��Ϊ������
				//ע�⣺1�����3�����������2�仯������
				//......2��z���������Ǻ�ɫ�ģ�����������1��
				//�Բ����£�a,b,c
				parent->color = BLACK;                           //a. �����ͼΪ��ɫ

				gparent->color = RED;                            //b. �游����Ϊ��ɫ

				root = rb_rotate_right(gparent, root);             //c. ���游���Ϊ֧������
			}
		}
		else    //if(parent == gparent->right)�����������Ϊ�游�����Һ���ʱ
		{
			uncle = gparent->left;                           //�游����������Ϊ������
			if (uncle && uncle->color == RED)                    //���1��z������y�Ǻ�ɫ��
			{
				uncle->color = BLACK;
				parent->color = BLACK;
				gparent->color = RED;
				node = gparent;

			}
			else                                    //���2��z������y�Ǻ�ɫ��
			{
				if (parent->left == node)                 //��zΪ����
				{
					root = rb_rotate_right(parent, root);        //����
					tmp = parent;
					parent = node;
					node = tmp;                      //������ɫ
				}
				//�������2�ı仯����Ϊ�����3.
				parent->color = BLACK;
				gparent->color = RED;
				root = rb_rotate_left(gparent, root);
			}
		}
	}

	root->color = BLACK;                             //����㣬������ô����������Ϊ��ɫ
	return root;                                    //���ظ����
}

//�����������
rb_node_t *rb_insert(key_t key, data_t data, rb_node_t *root)
{
	rb_node_t *parent = NULL, *node;
	parent = NULL;
	if ((node = rb_search_auxiliary(key, root, &parent)))             //����rb_search_auxiliary���ҵ�������ĵط�
	{
		return root;
	}

	node = rb_new_node(key, data);                        //������
	node->parent = parent;                             //��ʼ��
	node->left = node->right = NULL;
	node->color = RED;

	if (parent)                                 //parent�������Ҫ������ĵط�
	{
		if (parent->key > key)
		{
			parent->left = node;
		}
		else
		{
			parent->right = node;
		}
	}
	else
	{
		root = node;
	}
	return rb_insert_rebalance(node, root);                         //������󣬵���rb_search_rebalance�޸������
}

//�������4��ɾ�����
/* ɾ�����֮�����µ���������Ľṹ */
//x��ʾҪɾ���Ľ�㣬*other��w��ʾ�ֵܽ��
static rb_node_t *rb_erase_rebalance(rb_node_t *node, rb_node_t *parent, rb_node_t *root)
{
	rb_node_t *other, *o_left, *o_right;                             //x���ֵ�*other, �ֵ�����*o_left, �Һ���*o_right
	while ((!node || node->color == BLACK) && node != root)
	{
		if (parent->left == node)
		{
			other = parent->right;
			if (other->color == RED)                  //���1��x���ֵܽ��w�Ǻ�ɫ��
			{                           //�Բ�����
				other->color = BLACK;                //a. ���ֵܽ��Ⱦ�ɺ�ɫ

				parent->color = RED;             //b. �Ѹ����Ⱦ�ɺ�ɫ
				root = rb_rotate_left(parent, root);     //��p[x]��������
				other = parent->right;                //x�����ֵܽ��new w ����ת֮ǰw��ĳ�����ӡ���ʵ���������Ч��
			}
			//���2��x���ֵܽ��w�Ǻ�ɫ����w����������ȫΪ��ɫ
			if ((!other->left || other->left->color == BLACK) && (!other->right || other->right->color == BLACK))
			{
				//�Բߣ��ѵ�ǰ�����ֵܽ���г�ȡ��һ�غ�ɫ׷�ӵ�������ϣ��Ѹ���㵱���µĵ�ǰ���
				other->color = RED;              //a. �ֵܽ��Ϊ��ɫ
				node = parent;                  //b. �Ѹ���㵱���µĵ�ǰ���

				parent = node->parent;
			}
			else
			{
				//���3��x���ֵܽ��w�Ǻ�ɫ�ģ���w�������Ǻ�ɫ���Һ���Ϊ��ɫ
				if (!other->right || other->right->color == BLACK)
				{
					if ((o_left = other->left))           //w��������left[w]����ɫ����    
					{
						o_left->color = BLACK;           //w��������ɫ�ɺ�->��
					}
					other->color = RED;              //�ֵܽ��w�ɺ�->��
					root = rb_rotate_right(other, root);            //��w�����������Ӷ�������ʻָ� 
					other = parent->right;               //�仯�󣬸������Һ��ӣ���Ϊ�µ��ֵܽ��
				}

				//���4��x���ֵܽ��w�Ǻ�ɫ�ģ��ֵܽ����Һ����Ǻ�ɫ
				other->color = parent->color;                 //a. ���ֵܽ��Ⱦ�ɵ�ǰ��㸸������ɫ
				parent->color = BLACK;                   //b. �ѵ�ǰ���ĸ����Ⱦ�ɺ�ɫ
				if (other->right)
				{
					other->right->color = BLACK;              //c. �ֵܽ����Һ���Ⱦ�ɺ�ɫ
				}

				root = rb_rotate_left(parent, root);                    //d. �Ե�ǰ���ĸ����Ϊ֧���������
				node = root;                     //e. �ѵ�ǰ���x��Ϊ��root
				//��ʱ�㷨������������������ʵ�ת��ȷ

				break;
			}
		}
		else    //������������������ԭ��һ�¡�������
		{
			other = parent->left;
			if (other->color == RED)
			{
				other->color = BLACK;
				parent->color = RED;
				root = rb_rotate_right(parent, root);
				other = parent->left;
			}
			if ((!other->left || other->left->color == BLACK) && (!other->right || other->right->color == BLACK))
			{
				other->color = RED;
				node = parent;
				parent = node->parent;
			}
			else
			{
				if (!other->left || other->left->color == BLACK)
				{
					if ((o_right = other->right))
					{
						o_right->color = BLACK;
					}
					other->color = RED;
					root = rb_rotate_left(other, root);
					other = parent->left;
				}

				other->color = parent->color;
				parent->color = BLACK;
				if (other->left)
				{
					other->left->color = BLACK;
				}
				root = rb_rotate_right(parent, root);
				node = root;
				break;
			}
		}
	}

	if (node)
	{
		node->color = BLACK;
	}

	return root;
}

//�������ɾ�����
rb_node_t *rb_erase(key_t key, rb_node_t *root)
{
	rb_node_t *child, *parent, *old, *left, *node;
	color_t color;
	child = NULL;
	if (!(node = rb_search_auxiliary(key, root, NULL)))  //����rb_search_auxiliary����Ҫɾ���Ľ�� 
	{
		printf("key %d is not exist!\\n");
		return root;
	}
	old = node;
	//ɾ��������������ӽ��
	//���ԣ�����ѡ��ѵ�ǰɾ������������е����Ԫ�ػ������������е���СԪ�طŵ��ŵ���ɾ������λ��
	if (node->left && node->right)
	{
		node = node->right;
		while ((left = node->left) != NULL)               //�ҵ�����������СԪ�ؽ�㣬���浽node��
		{
			node = left;
		}
		child = node->right;                     //childΪnode���Ҷ��ӣ���node����ӵ��ֵܽ��              
		parent = node->parent;                       //parentΪnode���ĸ����
		color = node->color;                     //color�������node����color
		if (child)                           //node���Ҷ��Ӳ�Ϊ�գ�������游�����Ϊ�����
		{
			child->parent = parent;
		}
		if (parent)                          //���node���ĸ���㲻Ϊ�գ�����Ϊ�����ɣ�
		{
			if (parent->left == node)                 //���н����ƶ�������node�ĺ��ӽ��ҽӵ�node�ĸ������   
			{
				parent->left = child;
			}
			else
			{
				parent->right = child;
			}                             //ͬ��
		}
		else                                   //node�ĸ����Ϊ�����
		{
			root = child;
		}

		if (node->parent == old)                  //node�ĸ�������Ҫɾ���Ľ��
		{
			parent = node;
		}

		node->parent = old->parent;                   //���а���������СԪ�ؽ��ŵ���ǰɾ����㴦�Ĳ���
		node->color = old->color;
		node->right = old->right;
		node->left = old->left;

		if (old->parent)                  //�����ǰɾ�����ĸ������ڣ���old��Ϊ�����
		{
			if (old->parent->left == old)          //��ǰɾ��������丸��������
			{
				old->parent->left = node;
			}
			else
			{
				old->parent->right = node;        //��ǰɾ��������丸�����Һ���
			}
		}
		else
		{
			root = node;                    //��ǰɾ�����Ϊ�����
		}

		old->left->parent = node;
		if (old->right)
		{
			old->right->parent = node;
		}
	}//if(node->left && node->right) 

	else        //��else��䴦����ǣ�ɾ�����û�ж��ӡ�ֻ��һ������ʱ�����       
	{
		if (!node->left)          //ɾ�����
		{
			child = node->right;
		}
		else if (!node->right)
		{
			child = node->left;
		}
		parent = node->parent;
		color = node->color;

		if (child)
		{
			child->parent = parent;
		}
		if (parent)
		{
			if (parent->left == node)
			{
				parent->left = child;
			}
			else
			{
				parent->right = child;
			}
		}
		else
		{
			root = child;
		}
	}

	free(old);
	if (color == BLACK)
	{
		root = rb_erase_rebalance(child, parent, root); //����rb_erase_rebalance���ָ����������
	}

	return root;
}

int _tmain(int argc, _TCHAR* argv[])
{
	int i, count = 7;
	key_t key;
	rb_node_t *root = NULL;
	rb_node_t *node = NULL;
	srand(time(NULL));
	for (i = 1; i < count; ++i)
	{
		key = rand() % count;
		if ((root = rb_insert(key, i, root)))
		{
			printf("[i = %d] insert key %d success!\n", i, key);
		}
		else
		{
			printf("[i = %d] insert key %d error!\n", i, key);
			exit(-1);
		}
		if ((node = rb_search(key, root)))
		{
			printf("[i = %d] search key %d success!\n", i, key);
		}
		else
		{
			printf("[i = %d] search key %d error!\n", i, key);
			exit(-1);
		}

		if (!(i % 5))
		{
			if ((root = rb_erase(key, root)))
			{
				printf("[i = %d] erase key %d success\n", i, key);
			}
			else
			{
				printf("[i = %d] erase key %d error\n", i, key);
			}
		}
	}
	system("pause");
	return 0;
}

