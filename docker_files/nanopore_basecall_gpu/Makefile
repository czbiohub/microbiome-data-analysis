NAME=nanopore_basecall_demux_gpu
TAG:=$(shell /bin/date +%Y%m%d%H%M%S)
REGISTRY=brianyu2010/$(NAME)

all: build push

build:
	docker image build -t $(NAME):$(TAG) -t $(NAME):latest -f Dockerfile . # add --no-cache after build if needed
	docker image tag $(NAME):latest $(REGISTRY):latest
	docker image tag $(NAME):$(TAG) $(REGISTRY):$(TAG)

push:
	docker image push $(REGISTRY):$(TAG)
	docker image push $(REGISTRY):latest
	echo "brianyu2010/$(NAME):$(TAG)" > LATEST

run:
	docker container run --runtime=nvidia --rm $(NAME):latest
