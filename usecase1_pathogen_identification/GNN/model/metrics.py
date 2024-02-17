import torch
from torchmetrics import Metric, Accuracy, AUROC, AveragePrecision, F1Score, PrecisionRecallCurve, ROC

class ModelMetrics(Metric):
    def __init__(self, num_classes, **kwargs):
        super().__init__(**kwargs)

        self.num_classes = num_classes
        if self.num_classes < 2:
            raise ValueError("num_classes must be greater than or equal to 2")
        
        if self.num_classes == 2:
            self.accuracy = Accuracy(task="binary")
            self.auroc = AUROC(task="binary")
            self.average_precision = AveragePrecision(task="binary")
            self.f1_score = F1Score(task="binary")
            self.pr_curve = PrecisionRecallCurve(task="binary")
            self.roc = ROC(task="binary")
        else:
            self.accuracy = Accuracy(task="multiclass", num_classes=self.num_classes)
            self.auroc = AUROC(task="multiclass", num_classes=self.num_classes)
            self.average_precision = AveragePrecision(task="multiclass", num_classes=self.num_classes)
            self.f1_score = F1Score(task="multiclass", num_classes=self.num_classes)
            self.pr_curve = PrecisionRecallCurve(task="multiclass", num_classes=self.num_classes)
            self.roc = ROC(task="multiclass", num_classes=self.num_classes)

    def reset(self) -> None:
        self.accuracy.reset()
        self.auroc.reset()
        self.average_precision.reset()
        self.f1_score.reset()
        self.pr_curve.reset()
        self.roc.reset()

    def update(self, preds, target) -> None:
        # Ensure that preds and target are valid
        if not preds.dtype in (torch.float32, torch.float64):
            raise ValueError("preds must be a tensor of torch.float32 or torch.float64")
        
        if not target.dtype in (torch.int32, torch.int64):
            raise ValueError("target must be a tensor of torch.int32 or torch.int64")
        
        if self.num_classes == 2:
            if preds.shape != target.shape:
                raise ValueError("preds and target must have the same shape")
        else:
            if len(preds.shape) != 2 or preds.shape[1] != self.num_classes:
                raise ValueError("preds must be a tensor of shape (N, C) where C is the number of classes")
            if len(target.shape) != 1:
                raise ValueError("target must be a tensor of shape (N,)")

        # Update the metrics
        self.accuracy(preds, target)
        self.auroc(preds, target)
        self.average_precision(preds, target)
        self.f1_score(preds, target)
        self.pr_curve(preds, target)
        self.roc(preds, target)

    def compute(self):
        # Compute each metric
        accuracy = self.accuracy.compute()
        auroc = self.auroc.compute()
        average_precision = self.average_precision.compute()
        f1_score = self.f1_score.compute()
        pr_curve = self.pr_curve.compute()
        roc = self.roc.compute()

        return {'accuracy': accuracy.numpy(), 'auroc': auroc.numpy(), 'average_precision': average_precision.numpy(), 'f1_score': f1_score.numpy(), 'precision': pr_curve[0].numpy(), 'recall': pr_curve[1].numpy(), 'fpr': roc[0].numpy(), 'tpr': roc[1].numpy()}