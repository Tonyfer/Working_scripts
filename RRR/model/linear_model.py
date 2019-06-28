def mean_square_error(y_true, y_pred, sample_weights=None):
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    assert len(y_true) == len(y_pred)
        
    if type(sample_weights) == type(None):
        return(np.mean((y_true - y_pred)**2))
    else:
        sample_weights = np.array(sample_weights)
        assert len(sample_weights) == len(y_true)
        return(np.dot(sample_weights, (y_true - y_pred)**2))
    
loss_function = mean_square_error








class CustomLinearModel:
    """
    Linear model: Y = XB, fit by minimizing the provided loss_function
    with L2 regularization
    """
    def __init__(self, loss_function=mean_absolute_percentage_error, 
                 X=None, Y=None, sample_weights=None, reliability=None, beta_init=None, 
                 regularization=0.00012, theta = 1):
        self.reliability = reliability
        self.regularization = regularization
        self.theta = theta
        self.beta = None
        self.loss_function = loss_function
        self.sample_weights = sample_weights
        self.beta_init = beta_init
        
        self.X = X
        self.Y = Y
            
    
    def predict(self, X):
        prediction = np.matmul(X, self.beta)
        return(prediction)

    def model_error(self):
        error = self.loss_function(
            self.predict(self.X), self.Y, sample_weights=self.sample_weights
        )
        return(error)
    
    def l2_regularized_loss(self, beta):
        self.beta = beta
        return(self.model_error() + \
               sum(self.regularization*np.array(self.beta)**2))
               
    
    def secondary_regu_loss(self, beta):
        self.beta = beta
        return(self.l2_regularized_loss()-sum(self.theta*np.array(self.reliability)))
               
    
    
    def fit(self, maxiter=250):        
        # Initialize beta estimates (you may need to normalize
        # your data and choose smarter initialization values
        # depending on the shape of your loss function)
        if type(self.beta_init)==type(None):
            # set beta_init = 1 for every feature
            self.beta_init = np.array([1]*self.X.shape[1])
        else: 
            # Use provided initial values
            pass
            
        if self.beta!=None and all(self.beta_init == self.beta):
            print("Model already fit once; continuing fit with more itrations.")
        
        if type(reliability) == type(None):
            res = minimize(self.l2_regularized_loss, self.beta_init,
                           method='BFGS', options={'maxiter': 500})
        else:
            res = minimize(self.secondary_regu_loss, self.beta_init,
                           method='BFGS', options={'maxiter': 500})
        self.beta = res.x
        self.beta_init = self.beta

l2_mape_model = CustomLinearModel(
    loss_function=r,
    X=X, Y=Y, regularization=0.00012
)
l2_mape_model.fit()
l2_mape_model.beta