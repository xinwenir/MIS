import tensorflow as tf
import scipy.sparse as sp
import numpy as np


def loss_fn(A,x,b,d_err):
    # return tf.reduce_sum(tf.square(tf.sparse.sparse_dense_matmul(A, x) - b)) #仿真认为derr为1
  return  tf.reduce_sum(tf.square(tf.divide(tf.sparse.sparse_dense_matmul(A, x) - b,d_err))) 

# 定义罚函数
def penalty_function(x, lower_bound, upper_bound, penalty_factor):
    lower_penalty =  tf.maximum(lower_bound - x, 0)  # 下界惩罚项
    upper_penalty =  tf.maximum(x - upper_bound, 0)  # 上界惩罚项
    return penalty_factor * tf.reduce_sum(tf.square(lower_penalty + upper_penalty))

def Solver_equation(A,b,x0,bounds,d_err,refs):
    filename=r"E:\latex\my_paper\code\data\log"
    old_loss=None
    b = tf.constant([[i] for i in b],dtype=tf.float32)
    lower_bound=tf.constant([[i[0]] for i in bounds],dtype=tf.float32)
    upper_bound=tf.constant([[i[1]] for i in bounds],dtype=tf.float32)
    d_err=tf.constant([[i] for i in d_err],dtype=tf.float32)
    # lower_bound=tf.constant([[0] for i in bounds],dtype=tf.float32)
    # upper_bound=tf.constant([[2.65] for i in bounds],dtype=tf.float32)
    # 定义变量x0
    x = tf.Variable([[i] for i in x0.tolist()],dtype=tf.float32)
    refs = tf.Variable([[i] for i in refs],dtype=tf.float32)
    learning_rate=0.003
    beta=0.01
    penalty_factor=5
    filename+=str(learning_rate)+"_"+str(beta)+"_"+str(penalty_factor)

    # log_file=
    #定义学习率调度器
    initial_learning_rate = learning_rate
    decay_steps = 10
    decay_rate = 1
    learning_rate_fn = tf.keras.optimizers.schedules.ExponentialDecay(
        initial_learning_rate, decay_steps, decay_rate)
    
    
    # 定义损失函数和优化器
    optimizer = tf.keras.optimizers.SGD(name="tf",learning_rate=learning_rate_fn)
    # 迭代更新参数
    for i in range(0,500):
        # 计算损失函数和梯度
        with tf.GradientTape() as tape:
            loss1 = loss_fn(A,x,b,d_err)
            loss2 = penalty_function(x,lower_bound=lower_bound,upper_bound=upper_bound,penalty_factor=penalty_factor)
            loss3 = beta* tf.norm(tf.subtract(x,refs))**2
            loss=loss1+loss2+loss3
            gradients = tape.gradient(loss, x)
        
        # # #手动计算梯度和损失
        # loss,gradients=loss_grad_fn(A,x,b,bounds,0.001)
        
        # clipped_gradients = tf.clip_by_value(gradients, lower_bound, upper_bound)  # 对梯度进行投影
        # # 对变量进行投影到约束域内
        # x.assign(tf.clip_by_value(x, lower_bound, upper_bound))
        # 更新参数
        optimizer.apply_gradients([(gradients, x)])

        # 打印损失函数值
        if i % 10 == 0:
            # x.assign(tf.clip_by_value(x, lower_bound, upper_bound))
            if i%10==0:
                print("Iteration {},学习率：{}, misfit: {:.2e},罚函数: {:.2e},norm:{:.2e},loss: {:.2e}".format(i,optimizer.learning_rate.numpy() ,loss1.numpy(),loss2.numpy(),loss3.numpy(),loss.numpy()))
                if old_loss is not None and abs(old_loss-loss)<0.1:
                    break
    # x.assign(tf.clip_by_value(x, lower_bound, upper_bound))
    res=[]
    middle=x.numpy()
    for d in middle:
        res.append(d[0])
    
    return res

